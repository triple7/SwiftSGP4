//
//  File.swift
//  
//
//  Created by Yuma decaux on 25/6/2023.
//

import Foundation
import simd

extension SwiftSGP4 {

    func dot(_ x: [Double], _ y: [Double]) -> Double {
        return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]
    }

    func sgn(_ x: Double) -> Double {
        if x < 0.0 {
            return -1.0
        } else {
            return 1.0
        }
    }

    func mag(_ x: [Double]) -> Double {
        return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
    }

    func angle(_ vec1: [Double], _ vec2: [Double]) -> Double {
        let small = 0.00000001
        let undefined = 999999.1
        
        let magv1 = mag(vec1)
        let magv2 = mag(vec2)
        
        if magv1 * magv2 > small * small {
            var temp = dot(vec1, vec2) / (magv1 * magv2)
            if abs(temp) > 1.0 {
                temp = sgn(temp) * 1.0
            }
            return acos(temp)
        } else {
            return undefined
        }
    }

    
    func jday(year: Int, mon: Int, day: Int, hr: Int, minute: Int, sec: Double, timezone: Int, daylightsaving: Bool) -> Double {
        let a: Double = 367.0 * Double(year)
        let b: Double = floor((7.0 * (Double(year) + floor(Double(mon + 9) / 12.0))) * 0.25)
        let c: Double = floor(275.0 * Double(mon) / 9.0)
        let d: Double = Double(day) + 1721013.5

        let minutesInHour: Double = Double(minute) / 60.0
        let secondsInMinute: Double = sec / 60.0
        let daylightSummer:Double = (daylightsaving && summertime(year: year, month: mon, day: day, hour: hr, timezone: timezone)) ? 1.0 : 0.0
        let timeOffset: Double = Double(hr) - Double(timezone) - daylightSummer
        let e: Double = (secondsInMinute + minutesInHour + timeOffset) / 24.0
        
        let jd: Double = a - b + c + d + e
        return jd
    }


    func gstime(jdut1: Double) -> Double {
        let twopi = 2.0 * .pi
        let deg2rad = .pi / 180.0
        var temp, tut1: Double
        
        tut1 = (jdut1 - 2451545.0) / 36525.0
        temp = -6.2e-6 * tut1 * tut1 * tut1 + 0.093104 * tut1 * tut1 + (876600.0 * 3600 + 8640184.812866) * tut1 + 67310.54841 // sec
        temp = floatmod(a: temp * deg2rad / 240.0, b: twopi) // 360/86400 = 1/240, to deg, to rad

        // Check quadrants
        if temp < 0.0 {
            temp += twopi
        }

        return temp
    }

    func days2mdhms(year: Int, days: Double) -> (mon: Int, day: Int, hr: Int, minute: Int, sec: Double) {
        var mon = 0, day = 0, hr = 0, minute = 0, sec = 0.0
        let dayofyr = Int(floor(days))
        var lmonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        
        if (year % 4) == 0 {
            lmonth[1] = 29
        }
        
        var i = 0
        var inttemp = 0
        
        while (dayofyr > inttemp + lmonth[i]) && (i < 12) {
            inttemp += lmonth[i]
            i += 1
        }
        
        mon = i + 1
        day = dayofyr - inttemp
        
        var temp = (days - Double(dayofyr)) * 24.0
        hr = Int(floor(temp))
        
        temp = (temp - Double(hr)) * 60.0
        minute = Int(floor(temp))
        sec = (temp - Double(minute)) * 60.0
        
        return (mon, day, hr, minute, sec)
    }

    func invjday(jd: Double, timezone: Int, daylightsaving: Bool) -> (year: Int, mon: Int, day: Int, hr: Int, minute: Int, sec: Double) {
        var year = 0, mon = 0, day = 0, hr = 0, minute = 0, sec = 0.0
        var leapyrs: Int
        var days, tu, temp: Double
        
        let jdAdjusted = jd + Double(timezone) / 24.0
        
        // Find year and days of the year
        temp = jdAdjusted - 2415019.5
        tu = temp / 365.25
        year = 1900 + Int(floor(tu))
        leapyrs = Int(floor(Double(year - 1901) * 0.25))
        
        // Optional nudge by 8.64x10^-7 sec to get even outputs
        days = temp - (Double(year - 1900) * 365.0 + Double(leapyrs)) + 0.00000000001
        
        // Check for the case of the beginning of a year
        if days < 1.0 {
            year -= 1
            leapyrs = Int(floor(Double(year - 1901) * 0.25))
            days = temp - (Double(year - 1900) * 365.0 + Double(leapyrs))
        }
        
        // Find remaining data
        let result = days2mdhms(year: year, days: days)
        mon = result.mon
        day = result.day
        hr = result.hr
        minute = result.minute
        sec = result.sec
        
        if daylightsaving && summertime(year: year, month: mon, day: day, hour: hr, timezone: timezone) {
            let result2 = days2mdhms(year: year, days: days + 1.0/24.0)
            mon = result2.mon
            day = result2.day
            hr = result2.hr
            minute = result2.minute
            sec = result2.sec
        }
        
        sec -= 0.00000086400
        
        return (year, mon, day, hr, minute, sec)
    }

    func floatmod(a: Double, b: Double) -> Double {
        return a - b * floor(a / b)
    }

    
    func summertime(year: Int, month: Int, day: Int, hour: Int, timezone: Int) -> Bool {
        if month < 3 || month > 10 {
            return false // No daylight saving time in Jan, Feb, Nov, Dec
        }
        if month > 3 && month < 10 {
            return true // Daylight saving time in Apr, May, Jun, Jul, Aug, Sep
        }
        
        let leapYear = (5 * year / 4 + 4) % 7
        
        if (month == 3 && (hour + 24 * day) >= (1 + timezone + 24 * (31 - leapYear))) ||
           (month == 10 && (hour + 24 * day) < (1 + timezone + 24 * (31 - ((5 * year / 4 + 1) % 7)))) {
            return true
        } else {
            return false
        }
    }

    
    func teme2Ecef(rteme: simd_double3, jdut1: Double) -> simd_double3 {
        var recef = simd_double3()
        
        let gmst = gstime(jdut1: jdut1)
        
        let st = simd_double3x3([
            simd_double3(cos(gmst), -sin(gmst), 0.0),
            simd_double3(sin(gmst), cos(gmst), 0.0),
            simd_double3(0.0, 0.0, 1.0)
        ])
        
        let rpef = simd_mul(st, rteme)
        
        let pm = polarm(jdut1: jdut1)
        
        recef = simd_mul(pm, rpef)
        
        return recef
    }

    func polarm(jdut1: Double) -> simd_double3x3 {
        var pm = simd_double3x3()
        let pi = Double.pi
        let MJD: Double // Julian Date - 2,400,000.5 days
        let A: Double
        let C: Double
        let xp: Double // Polar motion coefficient in radians
        let yp: Double // Polar motion coefficient in radians
        
        // Predict polar motion coefficients using IERS Bulletin - A (Vol. XXVIII No. 030)
        MJD = jdut1 - 2400000.5
        A = 2 * pi * (MJD - 57226) / 365.25
        C = 2 * pi * (MJD - 57226) / 435
        
        xp = (0.1033 + 0.0494 * cos(A) + 0.0482 * sin(A) + 0.0297 * cos(C) + 0.0307 * sin(C)) * 4.84813681e-6
        yp = (0.3498 + 0.0441 * cos(A) - 0.0393 * sin(A) + 0.0307 * cos(C) - 0.0297 * sin(C)) * 4.84813681e-6
        
        pm[0][0] = cos(xp)
        pm[0][1] = 0.0
        pm[0][2] = -sin(xp)
        pm[1][0] = sin(xp) * sin(yp)
        pm[1][1] = cos(yp)
        pm[1][2] = cos(xp) * sin(yp)
        pm[2][0] = sin(xp) * cos(yp)
        pm[2][1] = -sin(yp)
        pm[2][2] = cos(xp) * cos(yp)
        
        return pm
    }

    func ijk2ll(r: SIMD3<Double>) -> SIMD3<Double> {
        var latlongh = SIMD3<Double>()
        let twopi = 2.0 * Double.pi
        let small = 0.00000001 // small value for tolerances
        let re = 6378.137 // radius of the Earth in km
        let eesqrd = 0.006694385000 // eccentricity of Earth squared
        
        let magr = simd.length(r)
        let temp = sqrt(r.x * r.x + r.y * r.y)
        
        if fabs(temp) < small {
            latlongh.y = simd.sign(r.z) * Double.pi * 0.5
        } else {
            latlongh.y = atan2(r.y, r.x)
        }
        
        if fabs(latlongh.y) >= Double.pi {
            if latlongh.y < 0.0 {
                latlongh.y += twopi
            } else {
                latlongh.y -= twopi
            }
        }
        
        latlongh.x = asin(r.z / magr)
        
        // Iterate to find geodetic latitude
        var i = 1
        var olddelta = latlongh.x + 10.0
        var sintemp: Double = 0
        var c: Double = 0
        
        while (fabs(olddelta - latlongh.x) >= small) && (i < 10) {
            olddelta = latlongh.x
            sintemp = sin(latlongh.x)
            c = re / sqrt(1.0 - eesqrd * sintemp * sintemp)
            latlongh.x = atan((r.z + c * eesqrd * sintemp) / temp)
            i += 1
        }
        
        if 0.5 * Double.pi - fabs(latlongh.x) > Double.pi / 180.0 {
            latlongh.z = (temp / cos(latlongh.x)) - c
        } else {
            latlongh.z = r.z / sin(latlongh.x) - c * (1.0 - eesqrd)
        }
        
        return latlongh
    }

    func site(latgd: Double, lon: Double, alt: Double) -> SIMD3<Double> {
        var rs = SIMD3<Double>()
        let re = 6378.137 // radius of the Earth in km
        let eesqrd = 0.006694385000 // eccentricity of Earth squared
        
        // Find rdel and rk components of the site vector
        let sinlat = sin(latgd)
        let cearth = re / sqrt(1.0 - (eesqrd * sinlat * sinlat))
        let rdel = (cearth + alt) * cos(latgd)
        let rk = ((1.0 - eesqrd) * cearth + alt) * sinlat
        
        // Find site position vector (ECEF)
        rs.x = rdel * cos(lon)
        rs.y = rdel * sin(lon)
        rs.z = rk
        
        return rs
    }

    public func rot3(invec: SIMD3<Double>, xval: Double) -> SIMD3<Double> {
        var outvec = SIMD3<Double>()
        let temp = invec.y
        let c = cos(xval)
        let s = sin(xval)
        
        outvec.y = c * invec.y - s * invec.x
        outvec.x = c * invec.x + s * temp
        outvec.z = invec.z
        
        return outvec
    }

    public func rot2(invec: SIMD3<Double>, xval: Double) -> SIMD3<Double> {
        var outvec = SIMD3<Double>()
        let temp = invec.z
        let c = cos(xval)
        let s = sin(xval)
        
        outvec.z = c * invec.z + s * invec.x
        outvec.x = c * invec.x - s * temp
        outvec.y = invec.y
        
        return outvec
    }

    
}
