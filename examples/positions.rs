//! Quick start: positions, an ayanamsha, houses, a star and Mercury's stations.
//! Run from the crate root with the ephemeris files in `data/`:
//! `cargo run --release --example positions`
use astro_calc_core::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let eph = Ephemeris::open_dir("data")?;
    let jd = JdUt::from_iso("1996-01-12T11:25:00")?;
    let geo = GeoPos::new(56.1, -3.6, 0.0)?;

    let bodies = parse_body_list("su,mo,me,ve,ma,ju,sa,ur,ne,pl,ra,ml")?;
    for p in eph.calc_positions(&bodies, jd, geo, CoordSystem::Ecliptic, Perspective::Topocentric)? {
        println!("{:3} {:9.4} {:+8.4}  {:+.4} deg/day", p.body.to_string(), p.lng, p.lat, p.lng_speed);
    }

    let lahiri = eph.ayanamsha(jd, Ayanamsha::Lahiri)?;
    println!("Lahiri ayanamsha {lahiri:.4}");

    let angles = angles(jd, geo);
    println!("Ascendant {:.4}, MC {:.4}", angles.ascendant, angles.mc);
    println!("Placidus cusps {:.2?}", houses(jd, geo, HouseSystem::Placidus).cusps);

    if let Some(spica) = star_by_name("Spica") {
        let p = eph.calc_position(spica.body(), jd, geo, CoordSystem::Ecliptic, Perspective::Geocentric)?;
        println!("Spica {:.4}", p.lng);
    }

    let mercury = eph.next_stations(Body::Mercury, jd, 4, 0)?;
    println!("Mercury retrograde: {}", mercury.is_retrograde());
    for s in &mercury.stations {
        println!("  {:11} {}", s.kind().key(), JdUt(s.jd()).to_iso()?);
    }
    Ok(())
}
