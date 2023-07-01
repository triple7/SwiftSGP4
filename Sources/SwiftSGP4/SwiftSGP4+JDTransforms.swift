//
//  File.swift
//  
//
//  Created by Yuma decaux on 22/6/2023.
//

import Foundation
import SGP4

let unix2000 = utc2000.timeIntervalSince1970 * 1000
let calendar = Calendar(identifier: .gregorian)
let utc2000Components = DateComponents(year: 2000, month: 1, day: 1)
let jd2000:Double = 2451545.0
let utc2000 = calendar.date(from: utc2000Components)!

extension SwiftSGP4 {
    
    
    func getJulianFromUnix(_ unixSecs: Double) -> Double {
        return (unixSecs / 86400.0) + 2440587.5
    }

    func getUnixFromJulian(_ julian: Double) -> UInt64 {
        return UInt64((julian - 2440587.5) * 86400.0 + 0.5)
    }

    func date2jday(_ dt: Date) -> Double {
        let dateUnix = dt.timeIntervalSince1970 * 1000
        let deltaMSec = (dateUnix - unix2000)
        let jd = jd2000 + deltaMSec / 86400000.0
        return jd
    }

    func jday2date(_ jd: Double) -> Date {
        let jdDeltaMSec = (jd - jd2000) * 86400000.0
        let date = addMilliseconds(utc2000, milliseconds: jdDeltaMSec)
        return date
    }

    func addDays(_ dt: Date, days: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.day = days
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    func addMinutes(_ dt: Date, minutes: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.minute = minutes
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    func addSeconds(_ dt: Date, seconds: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.second = seconds
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    func addMilliseconds(_ dt: Date, milliseconds: Double) -> Date {
        let updatedDate = dt.addingTimeInterval(milliseconds / 1000)
        return updatedDate
    }
}
