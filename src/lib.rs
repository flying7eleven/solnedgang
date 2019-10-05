use chrono::{DateTime, Datelike, Utc};
use std::fmt;

/// Representation of a location on earth.
pub struct Coordinates {
    latitude: f64,
    longitude: f64,
}

impl Coordinates {
    /// Create a new `Location` instance from a pair of latitude and longitude in decimal representation.
    ///
    /// # Examples
    ///
    /// To create a new instance of this class and print the coordinates to the console:
    /// ```
    /// use solnedgang::Coordinates;
    ///
    /// let duesseldorf = Coordinates::new(51.21794, 6.76165);
    /// println!("GPS coordinates for Düsseldorf, Germany are: {}", duesseldorf);
    /// ```
    pub fn new(latitude: f64, longitude: f64) -> Coordinates {
        Coordinates {
            latitude,
            longitude,
        }
    }
}

impl fmt::Display for Coordinates {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "({}, {})", self.latitude, self.longitude)
    }
}

/// Calculator for determine the sunrise and sunset times.
pub struct SunriseSunsetCalculator {
    location: Coordinates,
}

impl From<Coordinates> for SunriseSunsetCalculator {
    /// Create a new instance of the `SunriseSunsetCalculator` class based on the location supplied.
    ///
    /// # Examples
    ///
    /// Create the calculator for Düsseldorf, Germany:
    /// ```
    /// use solnedgang::{SunriseSunsetCalculator, Coordinates};
    ///
    /// let duesseldorf = Coordinates::new(51.21794, 6.76165);
    /// let calculator = SunriseSunsetCalculator::from(duesseldorf);
    /// ```
    fn from(location: Coordinates) -> Self {
        SunriseSunsetCalculator { location }
    }
}

impl SunriseSunsetCalculator {
    fn calc_mean_obliquity_of_ecliptic(&self, t: f64) -> f64 {
        let seconds = 21.448 - t * (46.815_0 + t * (0.000_59 - t * (0.001_813)));
        23.0 + (26.0 + (seconds / 60.0)) / 60.0
    }

    fn calc_geom_mean_long_sun(&self, t: f64) -> f64 {
        //
        let mut l = 280.466_46 + t * (36_000.769_83 + 0.000_303_2 * t);

        //
        loop {
            if l > 360.0 {
                l -= 360.0;
            } else {
                break;
            }
        }

        //
        loop {
            if l < 0.0 {
                l += 360.0;
            } else {
                break;
            }
        }

        //
        l
    }

    fn calc_obliquity_correction(&self, t: f64) -> f64 {
        let e0 = self.calc_mean_obliquity_of_ecliptic(t);
        let omega = 125.04 - 1934.136 * t;
        e0 + 0.002_56 * omega.to_radians().cos()
    }

    fn calc_eccentricity_earth_orbit(&self, t: f64) -> f64 {
        0.016_708_634 - t * (0.000_042_037 + 0.000_000_126_7 * t)
    }

    fn calc_geom_mean_anomaly_sun(&self, t: f64) -> f64 {
        357.529_11 + t * (35_999.050_29 - 0.000_153_7 * t)
    }

    fn calc_equation_of_time(&self, t: f64) -> f64 {
        let epsilon = self.calc_obliquity_correction(t);
        let l0 = self.calc_geom_mean_long_sun(t);
        let e = self.calc_eccentricity_earth_orbit(t);
        let m = self.calc_geom_mean_anomaly_sun(t);
        let mut y = (epsilon.to_radians() / 2.0).tan();

        y *= y;

        let sin2l0 = (2.0 * l0.to_radians()).sin();
        let sinm = m.to_radians().sin();
        let cos2l0 = (2.0 * l0.to_radians()).cos();
        let sin4l0 = (4.0 * l0.to_radians()).sin();
        let sin2m = (2.0 * m.to_radians()).sin();
        let etime = y * sin2l0 - 2.0 * e * sinm + 4.0 * e * y * sinm * cos2l0
            - 0.5 * y * y * sin4l0
            - 1.25 * e * e * sin2m;

        etime.to_degrees() * 4.0 // in minutes of time
    }

    fn calc_time_julian_cent(&self, jd: f64) -> f64 {
        (jd - 2_451_545.0) / 36_525.0
    }

    fn calc_sun_eq_of_center(&self, t: f64) -> f64 {
        let m = self.calc_geom_mean_anomaly_sun(t);
        let mrad = m.to_radians();
        let sinm = mrad.sin();
        let sin2m = (mrad + mrad).sin();
        let sin3m = (mrad + mrad + mrad).sin();

        sinm * (1.914_602 - t * (0.004_817 + 0.000_014 * t))
            + sin2m * (0.019_993 - 0.000_101 * t)
            + sin3m * 0.000_289
    }

    fn calc_sun_true_long(&self, t: f64) -> f64 {
        let l0 = self.calc_geom_mean_long_sun(t);
        let c = self.calc_sun_eq_of_center(t);

        l0 + c // in degrees
    }

    fn calc_sun_apparent_long(&self, t: f64) -> f64 {
        let o = self.calc_sun_true_long(t);
        let omega = 125.04 - 1934.136 * t;
        o - 0.005_69 - 0.004_78 * omega.to_radians().sin() // in degrees
    }

    fn calc_sun_declination(&self, t: f64) -> f64 {
        let e = self.calc_obliquity_correction(t);
        let lambda = self.calc_sun_apparent_long(t);

        let sint = e.to_radians().sin() * lambda.to_radians().sin();
        sint.asin().to_degrees() // in degrees
    }

    fn calc_hour_angle_sunrise(&self, lat: f64, solar_dec: f64) -> f64 {
        let lat_rad = lat.to_radians();
        let sd_rad = solar_dec.to_radians();
        let z: f64 = 90.833;
        (z.to_radians().cos() / (lat_rad.cos() * sd_rad.cos()) - lat_rad.tan() * sd_rad.tan())
            .acos() // in radians
    }

    fn calc_hour_angle_sunset(&self, lat: f64, solar_dec: f64) -> f64 {
        let lat_rad = lat.to_radians();
        let sd_rad = solar_dec.to_radians();
        let z: f64 = 90.833;
        -(z.to_radians().cos() / (lat_rad.cos() * sd_rad.cos()) - lat_rad.tan() * sd_rad.tan())
            .acos() // in radians
    }

    fn calc_jd(&self, year_in: i32, month_in: u32, day: u32) -> f64 {
        let mut year = year_in;
        let mut month = month_in;

        if month <= 2 {
            year -= 1;
            month += 12;
        }

        let a = (f64::from(year) / 100.0).floor();
        let b = 2.0 - a + (a / 4.0).floor();

        (365.25 * (f64::from(year) + 4716.0)).floor()
            + (30.600_1 * (f64::from(month) + 1.0)).floor()
            + f64::from(day)
            + b
            - 1_524.5
    }

    fn calc_jdfrom_julian_cent(&self, t: f64) -> f64 {
        t * 36_525.0 + 2_451_545.0
    }

    /// Calculates the the time for the sunrise in minutes from midnight (UTC).
    ///
    /// # Arguments
    /// * `date` - The date for which the sunrise time should be calculated for.
    pub fn calc_sunrise_utc(&self, date: DateTime<Utc>) -> f64 {
        let julian_date = self.calc_jd(date.year(), date.month(), date.day());
        let t = self.calc_time_julian_cent(julian_date);
        // first pass to approximate sunrise
        let mut eq_time = self.calc_equation_of_time(t);
        let mut solar_dec = self.calc_sun_declination(t);
        let mut hour_angle = self.calc_hour_angle_sunrise(self.location.latitude, solar_dec);
        let mut delta = self.location.longitude + hour_angle.to_degrees();
        let mut time_diff = 4.0 * delta; // in minutes of time
        let time_utc = 720.0 - time_diff - eq_time; // in minutes
        let newt = self.calc_time_julian_cent(self.calc_jdfrom_julian_cent(t) + time_utc / 1_440.0);

        eq_time = self.calc_equation_of_time(newt);
        solar_dec = self.calc_sun_declination(newt);

        hour_angle = self.calc_hour_angle_sunrise(self.location.latitude, solar_dec);
        delta = self.location.longitude + hour_angle.to_degrees();
        time_diff = 4.0 * delta;

        720.0 - time_diff - eq_time // return time in minutes from midnight
    }

    /// Calculates the the time for the sunset in minutes from midnight (UTC).
    ///
    /// # Arguments
    /// * `date` - The date for which the sunset time should be calculated for.
    pub fn calc_sunset_utc(&self, date: DateTime<Utc>) -> f64 {
        let julian_date = self.calc_jd(date.year(), date.month(), date.day());
        let t = self.calc_time_julian_cent(julian_date);
        // first pass to approximate sunset
        let mut eq_time = self.calc_equation_of_time(t);
        let mut solar_dec = self.calc_sun_declination(t);
        let mut hour_angle = self.calc_hour_angle_sunset(self.location.latitude, solar_dec);
        let mut delta = self.location.longitude + hour_angle.to_degrees();
        let mut time_diff = 4.0 * delta; // in minutes of time
        let mut time_utc = 720.0 - time_diff - eq_time; // in minutes
        let newt = self.calc_time_julian_cent(self.calc_jdfrom_julian_cent(t) + time_utc / 1_440.0);

        eq_time = self.calc_equation_of_time(newt);
        solar_dec = self.calc_sun_declination(newt);

        hour_angle = self.calc_hour_angle_sunset(self.location.latitude, solar_dec);
        delta = self.location.longitude + hour_angle.to_degrees();
        time_diff = 4.0 * delta;
        time_utc = 720.0 - time_diff - eq_time; // in minutes

        time_utc
        // time_utc + (60 * tzOffset) // return time in minutes from midnight
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;

    #[test]
    fn calculating_sunrise_in_utc_seconds_works_for_duesseldorf_germany() {
        let location = Coordinates::new(51.21875, 6.76341);
        let calculator = SunriseSunsetCalculator { location };

        let test_date = Utc.ymd(2019, 10, 4).and_hms(0, 0, 0);

        assert_eq!(
            (calculator.calc_sunrise_utc(test_date) - 337.668) > 0.0,
            true
        );
        assert_eq!(
            (calculator.calc_sunrise_utc(test_date) - 337.668) < 0.0001,
            true
        );
    }

    #[test]
    fn calculating_sunset_in_utc_seconds_works_for_duesseldorf_germany() {
        let location = Coordinates::new(51.21875, 6.76341);
        let calculator = SunriseSunsetCalculator { location };

        let test_date = Utc.ymd(2019, 10, 4).and_hms(0, 0, 0);

        assert_eq!(
            (calculator.calc_sunset_utc(test_date) - 1024.9458) > 0.0,
            true
        );
        assert_eq!(
            (calculator.calc_sunset_utc(test_date) - 1024.9458) < 0.0001,
            true
        );
    }
}
