//! Downloads and probes. Everything is written to `<name>.part`, checked (SPK index, the
//! bodies expected, SHA-256 when one is configured) and only then renamed into place, so a
//! failed run never leaves a broken `.bsp` for `Ephemeris::open_dir`.

use std::collections::HashMap;
use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Duration;

use anyhow::{Context, Result, anyhow, bail};
use astro_calc_core::{DateRange, Ephemeris};
use base64::Engine;
use indicatif::{ProgressBar, ProgressStyle};
use sha2::{Digest, Sha256};
use ureq::Agent;

pub fn agent() -> Agent {
    Agent::config_builder()
        .user_agent(concat!("astro-calc-build/", env!("CARGO_PKG_VERSION")))
        .timeout_connect(Some(Duration::from_secs(30)))
        // Horizons takes a while to generate long spans.
        .timeout_recv_response(Some(Duration::from_secs(600)))
        .http_status_as_error(false)
        .build()
        .into()
}

pub fn part_path(dest: &Path) -> PathBuf {
    let mut name = dest.file_name().unwrap_or_default().to_os_string();
    name.push(".part");
    dest.with_file_name(name)
}

fn file_name(p: &Path) -> String {
    p.file_name().unwrap_or_default().to_string_lossy().into_owned()
}

/// Downloads `url` to `dest`'s part file, resuming a previous partial download.
pub fn download(agent: &Agent, url: &str, dest: &Path, size_hint: Option<u64>) -> Result<PathBuf> {
    let part = part_path(dest);
    let label = file_name(dest);
    let have = fs::metadata(&part).map(|m| m.len()).unwrap_or(0);
    let mut req = agent.get(url);
    if have > 0 {
        req = req.header("Range", format!("bytes={have}-"));
    }
    let mut resp = req.call().with_context(|| format!("cannot reach {url}"))?;
    let status = resp.status().as_u16();
    let header = |name: &str| resp.headers().get(name).and_then(|v| v.to_str().ok()).map(str::to_owned);
    let (start, total) = match status {
        // Content-Range: bytes 1000-1999/2000
        206 => (have, header("content-range").and_then(|r| r.rsplit('/').next()?.parse().ok()).or(size_hint)),
        200 => (0, header("content-length").and_then(|n| n.parse().ok()).or(size_hint)),
        // The part file is already complete.
        416 if have > 0 => return Ok(part),
        s => bail!("{url}: HTTP {s}"),
    };
    if have > 0 && start == have {
        eprintln!("  resuming {label} at {}", crate::catalog::fmt_bytes(have as f64));
    }
    let out = if start > 0 { OpenOptions::new().append(true).open(&part)? } else { File::create(&part)? };
    let pb = match total {
        Some(t) => ProgressBar::new(t).with_style(
            ProgressStyle::with_template("  {msg:16} [{bar:30}] {bytes}/{total_bytes} {binary_bytes_per_sec} eta {eta}")
                .expect("template")
                .progress_chars("=> "),
        ),
        None => ProgressBar::new_spinner()
            .with_style(ProgressStyle::with_template("  {msg:16} {spinner} {bytes} {binary_bytes_per_sec}").expect("template")),
    };
    pb.set_message(label.clone());
    pb.set_position(start);
    let mut reader = pb.wrap_read(resp.body_mut().as_reader());
    let mut writer = BufWriter::with_capacity(1 << 20, out);
    let copied = std::io::copy(&mut reader, &mut writer);
    writer.flush()?;
    pb.finish_and_clear();
    copied.with_context(|| format!("download of {label} interrupted; run the same command again to resume"))?;
    let got = fs::metadata(&part)?.len();
    if let Some(t) = total
        && got != t
    {
        bail!("{label} is incomplete ({got} of {t} bytes); run the same command again to resume");
    }
    Ok(part)
}

/// Asks Horizons to generate an SPK file for `command` over `from`..`to` (`YYYY-MM-DD`),
/// writing it to `dest`'s part file.
pub fn horizons(agent: &Agent, api: &str, command: &str, from: &str, to: &str, dest: &Path) -> Result<PathBuf> {
    let pb = ProgressBar::new_spinner().with_style(ProgressStyle::with_template("  {spinner} {msg}").expect("template"));
    pb.set_message(format!("{}: JPL Horizons is generating {from} to {to}…", file_name(dest)));
    pb.enable_steady_tick(Duration::from_millis(120));
    let result = (|| -> Result<Vec<u8>> {
        let mut resp = agent
            .get(api)
            .query("format", "json")
            .query("COMMAND", format!("'{command}'"))
            .query("OBJ_DATA", "NO")
            .query("MAKE_EPHEM", "YES")
            .query("EPHEM_TYPE", "SPK")
            .query("START_TIME", format!("'{from}'"))
            .query("STOP_TIME", format!("'{to}'"))
            .call()
            .with_context(|| format!("cannot reach {api}"))?;
        let status = resp.status().as_u16();
        let mut raw = Vec::new();
        resp.body_mut().as_reader().read_to_end(&mut raw)?;
        let json: serde_json::Value = serde_json::from_slice(&raw)
            .map_err(|_| anyhow!("Horizons returned HTTP {status} and no JSON: {}", String::from_utf8_lossy(&raw[..raw.len().min(200)])))?;
        let Some(spk) = json.get("spk").and_then(|s| s.as_str()) else {
            let why = json.get("error").or_else(|| json.get("result")).and_then(|s| s.as_str()).unwrap_or("no file returned");
            bail!("Horizons refused it: {}", why.trim().lines().take(4).collect::<Vec<_>>().join(" "));
        };
        let clean: String = spk.chars().filter(|c| !c.is_ascii_whitespace()).collect();
        base64::engine::general_purpose::STANDARD.decode(clean).context("Horizons sent invalid base64")
    })();
    pb.finish_and_clear();
    let part = part_path(dest);
    fs::write(&part, result?)?;
    Ok(part)
}

/// What a downloaded file must hold.
pub enum Expect {
    /// The Sun, Moon and Earth.
    Planetary,
    /// This NAIF id.
    Target(i32),
}

/// Checks a part file and moves it into place. A part that fails is deleted, so the next
/// run starts afresh rather than resuming bad data.
pub fn install_part(part: &Path, dest: &Path, expect: Expect, sha256: Option<&str>) -> Result<u64> {
    let check = || -> Result<()> {
        if let Some(want) = sha256 {
            let got = sha256_file(part)?;
            if got != want {
                bail!("SHA-256 mismatch: expected {want}, got {got}");
            }
        }
        let eph = Ephemeris::open(part).context("not a usable SPK file")?;
        match expect {
            Expect::Planetary if !eph.files().first().is_some_and(|f| f.is_planetary()) => {
                bail!("does not hold the Sun, Moon and Earth")
            }
            Expect::Target(id) if !eph.targets().contains(&id) => bail!("does not hold NAIF id {id}"),
            _ => Ok(()),
        }
    };
    if let Err(e) = check() {
        let _ = fs::remove_file(part);
        return Err(e.context(format!("{} rejected", file_name(dest))));
    }
    fs::rename(part, dest)?;
    Ok(fs::metadata(dest)?.len())
}

pub fn sha256_file(path: &Path) -> Result<String> {
    let mut f = File::open(path)?;
    let mut hasher = Sha256::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finalize().iter().map(|b| format!("{b:02x}")).collect())
}

fn get_range(agent: &Agent, url: &str, from: u64, len: u64) -> Result<(Vec<u8>, Option<u64>)> {
    let mut resp = agent.get(url).header("Range", format!("bytes={from}-{}", from + len - 1)).call()?;
    if resp.status().as_u16() != 206 {
        bail!("{url}: the server does not serve byte ranges (HTTP {})", resp.status().as_u16());
    }
    let total = resp
        .headers()
        .get("content-range")
        .and_then(|v| v.to_str().ok())
        .and_then(|r| r.rsplit('/').next()?.parse().ok());
    let mut buf = Vec::new();
    resp.body_mut().as_reader().take(len).read_to_end(&mut buf)?;
    Ok((buf, total))
}

/// The size of a remote SPK file and the span shared by every body in it, read from its
/// index with a few small range requests, without downloading the file.
pub fn probe(agent: &Agent, url: &str) -> Result<(Option<u64>, Option<DateRange>)> {
    const REC: u64 = 1024;
    let (head, total) = match get_range(agent, url, 0, REC) {
        Ok(r) => r,
        // No byte ranges: the size, at least, from a HEAD request.
        Err(e) => {
            let resp = agent.head(url).call()?;
            let size = resp.headers().get("content-length").and_then(|v| v.to_str().ok()?.parse().ok());
            return if resp.status().is_success() && size.is_some() { Ok((size, None)) } else { Err(e) };
        }
    };
    if head.len() < REC as usize || &head[0..7] != b"DAF/SPK" || &head[88..96] != b"LTL-IEEE" {
        return Ok((total, None));
    }
    let i32_at = |b: &[u8], o: usize| i32::from_le_bytes(b[o..o + 4].try_into().expect("4 bytes"));
    let f64_at = |b: &[u8], o: usize| f64::from_le_bytes(b[o..o + 8].try_into().expect("8 bytes"));
    if (i32_at(&head, 8), i32_at(&head, 12)) != (2, 6) {
        return Ok((total, None));
    }
    let mut spans: HashMap<i32, (f64, f64)> = HashMap::new();
    let mut record = i32_at(&head, 76) as u64;
    for _ in 0..200 {
        if record == 0 {
            break;
        }
        let (rec, _) = get_range(agent, url, (record - 1) * REC, REC)?;
        if rec.len() < REC as usize {
            break;
        }
        let next = f64_at(&rec, 0) as u64;
        let nsum = (f64_at(&rec, 16) as usize).min(25);
        for i in 0..nsum {
            let w = (3 + i * 5) * 8;
            let (start, end, target) = (f64_at(&rec, w), f64_at(&rec, w + 8), i32_at(&rec, w + 16));
            let s = spans.entry(target).or_insert((f64::INFINITY, f64::NEG_INFINITY));
            *s = (s.0.min(start), s.1.max(end));
        }
        record = next;
    }
    if spans.is_empty() {
        return Ok((total, None));
    }
    const J2000: f64 = 2_451_545.0;
    let start = spans.values().map(|s| s.0).fold(f64::NEG_INFINITY, f64::max);
    let end = spans.values().map(|s| s.1).fold(f64::INFINITY, f64::min);
    Ok((total, Some(DateRange { start: J2000 + start / 86_400.0, end: J2000 + end / 86_400.0 })))
}
