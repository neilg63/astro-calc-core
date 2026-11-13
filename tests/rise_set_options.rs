//! Every rise/set option: disc alignment, refraction, twilight, horizon height, dip,
//! geocentric positions and projection onto the ecliptic.

use std::path::PathBuf;

use astro_calc_core::*;

fn eph() -> Option<Ephemeris> {
    let bsp = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/de440s.bsp");
    bsp.exists().then(|| Ephemeris::open(bsp).unwrap())
}

const SECOND: f64 = 1.0 / 86_400.0;

#[test]
fn every_option_moves_the_event_as_it_should() {
    let Some(eph) = eph() else { return };
    let london = GeoPos::new(51.5074, -0.1278, 0.0).unwrap();
    let midnight = JdUt::from_iso("2024-03-20T00:00:00").unwrap().value();
    let rise = |body: Body, geo: GeoPos, opts: RiseSetOptions| eph.next_transit(body, midnight, geo, TransitEvent::Rise, opts, 1.5).unwrap().unwrap().value();
    let set = |body: Body, geo: GeoPos, opts: RiseSetOptions| eph.next_transit(body, midnight, geo, TransitEvent::Set, opts, 1.5).unwrap().unwrap().value();
    let sun = |opts: RiseSetOptions| rise(Body::Sun, london, opts);
    let alt = |t: f64| eph.altitude(Body::Sun, t, london).unwrap();
    let base = RiseSetOptions::default();

    // The almanac sunrise for London on 20 March 2024 is 06:02 UTC.
    let almanac = sun(RiseSetOptions::upper_limb());
    assert_eq!(&JdUt(almanac).to_iso().unwrap()[..15], "2024-03-20T06:0", "{}", JdUt(almanac).to_iso().unwrap());
    let minutes = (almanac - midnight) * 1440.0 - 360.0;
    assert!((1.0..4.0).contains(&minutes), "sunrise {minutes} minutes after 06:00");

    // All twelve combinations of limb, refraction and fixed radius: the upper limb rises
    // first and the lower last; refraction brings every one of them forward.
    use DiscAlignment::*;
    for fixed in [false, true] {
        for refraction in [Refraction::Off, Refraction::Standard] {
            let t = |disc| sun(base.with_disc(disc).with_refraction(refraction).with_fixed_disc(fixed));
            assert!(t(UpperLimb) < t(Center) && t(Center) < t(LowerLimb));
            assert!((t(Center) - t(UpperLimb) - 100.0 * SECOND).abs() < 30.0 * SECOND, "limb to centre takes about 100 s");
        }
        for disc in [UpperLimb, Center, LowerLimb] {
            let t = |refraction| sun(base.with_disc(disc).with_refraction(refraction).with_fixed_disc(fixed));
            assert!(t(Refraction::Standard) < t(Refraction::Off) - 120.0 * SECOND);
        }
    }
    // The centre's geometric altitude at each event is the target the options define.
    assert!(alt(sun(base.with_refraction(Refraction::Off))).abs() < 1e-5);
    assert!((alt(sun(base)) + rise_set::HORIZON_REFRACTION_ARCMIN / 60.0).abs() < 1e-5);
    let upper_no_refraction = alt(sun(base.with_disc(UpperLimb).with_refraction(Refraction::Off)));
    assert!((upper_no_refraction + 0.267).abs() < 0.005, "upper limb, no refraction: {upper_no_refraction}");

    // Twilights begin in order before sunrise and end in order after sunset.
    let twilight = |h| RiseSetOptions::twilight(h);
    let (astro, nautical, civil) = (sun(twilight(Horizon::AstronomicalTwilight)), sun(twilight(Horizon::NauticalTwilight)), sun(twilight(Horizon::CivilTwilight)));
    assert!(astro < nautical && nautical < civil && civil < almanac);
    for (t, depth) in [(civil, -6.0), (nautical, -12.0), (astro, -18.0)] {
        assert!((alt(t) - depth).abs() < 1e-5, "altitude {} for {depth}", alt(t));
    }
    let dusk = |h| set(Body::Sun, london, twilight(h));
    let sunset = set(Body::Sun, london, RiseSetOptions::upper_limb());
    assert!(sunset < dusk(Horizon::CivilTwilight) && dusk(Horizon::CivilTwilight) < dusk(Horizon::NauticalTwilight));
    // Civil twilight at this latitude near the equinox lasts a little over half an hour.
    assert!(((almanac - civil) * 1440.0 - 33.0).abs() < 4.0, "civil dawn {} min before sunrise", (almanac - civil) * 1440.0);
    // The disc and refraction settings make no difference to a twilight.
    assert_eq!(civil, sun(twilight(Horizon::CivilTwilight).with_disc(LowerLimb).with_refraction(Refraction::Off).with_dip(true)));

    // A ridge two degrees high delays sunrise; the centre is then two degrees up.
    let ridge = sun(base.with_horizon(Horizon::Altitude(2.0)).with_refraction(Refraction::Off));
    assert!(ridge > sun(base.with_refraction(Refraction::Off)) + 600.0 * SECOND);
    assert!((alt(ridge) - 2.0).abs() < 1e-5);
    // Refraction is weaker higher up: about 18′ at 2° against 37′ on the horizon.
    let refracted = alt(sun(base.with_horizon(Horizon::Altitude(2.0))));
    assert!((2.0 - refracted) * 60.0 > 15.0 && (2.0 - refracted) * 60.0 < 22.0, "refraction at 2°: {}′", (2.0 - refracted) * 60.0);
    assert_eq!(sun(base.with_horizon(Horizon::Altitude(0.0))), sun(base));

    // Air: the standard atmosphere is 1010 hPa at 10 °C; thin or warm air refracts less.
    let air = |pressure, temperature| sun(base.with_refraction(Refraction::Atmosphere { pressure, temperature }));
    assert!((air(1010.0, 10.0) - sun(base)).abs() < 0.01 * SECOND);
    assert!(air(700.0, 10.0) > sun(base) + 30.0 * SECOND && air(1010.0, 35.0) > sun(base) + 5.0 * SECOND);
    assert!(air(1040.0, -20.0) < sun(base));

    // Dip: from 100 m up the horizon is 17.6′ lower, so the Sun rises about two minutes sooner.
    let tower = GeoPos::new(51.5074, -0.1278, 100.0).unwrap();
    let (flat, dipped) = (rise(Body::Sun, tower, base), rise(Body::Sun, tower, base.with_dip(true)));
    assert!(((flat - dipped) / SECOND - 110.0).abs() < 40.0, "dip advanced sunrise by {} s", (flat - dipped) / SECOND);
    let dip = eph.altitude(Body::Sun, flat, tower).unwrap() - eph.altitude(Body::Sun, dipped, tower).unwrap();
    assert!((dip * 60.0 - 17.6).abs() < 0.01, "dip {}′", dip * 60.0);
    assert_eq!(rise(Body::Sun, london, base.with_dip(true)), rise(Body::Sun, london, base), "no height, no dip");

    // Parallax matters for the Moon (about a degree) and hardly at all for the Sun.
    let moon_shift = (rise(Body::Moon, london, base.with_geocentric(true)) - rise(Body::Moon, london, base)) / SECOND;
    assert!(moon_shift < -120.0, "geocentric moonrise differs by {moon_shift} s");
    assert!(((sun(base.with_geocentric(true)) - sun(base)) / SECOND).abs() < 2.0);

    // Hindu rising: centre, no refraction, geocentric, on the ecliptic. The Sun is on the
    // ecliptic anyway; the Moon is not, so dropping its latitude moves the time.
    let hindu = RiseSetOptions::hindu_rising();
    assert_eq!((hindu.disc, hindu.refraction, hindu.geocentric, hindu.ignore_latitude), (Center, Refraction::Off, true, true));
    assert!(((sun(hindu) - sun(base.with_refraction(Refraction::Off).with_geocentric(true))) / SECOND).abs() < 1.0);
    let moon_lat = (rise(Body::Moon, london, hindu) - rise(Body::Moon, london, hindu.with_ignore_latitude(false))) / SECOND;
    assert!(moon_lat.abs() > 60.0, "ignoring the Moon's latitude moved moonrise by {moon_lat} s");
}
