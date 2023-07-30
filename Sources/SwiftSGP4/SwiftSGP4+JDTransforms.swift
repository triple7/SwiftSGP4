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
    
func dateString2Date( _ dateString: String)->Date {
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        let date = dateFormat.date(from: dateString)!
        return date
    }
    
    func timestampToJD( _ date: Date)->Double {
        let calendar = Calendar.current
        let  components  =  calendar.dateComponents([.year, .month, .day, .hour, .minute, .second],  from:  date)
        let year = components.year!.int32()
        let month = components.month!.int32()
        let day = components.day!.int32()
        let hour = components.hour!.int32()
        let minutes = components.minute!.int32()
         let seconds = Double(components.second!)
        var jd:Double = 0.0
        var jdFrac:Double = 0.0
        
        jday(year, month, day, hour, minutes, seconds, &jd, &jdFrac)
return jd + jdFrac
    }

extension SwiftSGP4 {
    
    public func framesSinceEpoch( _ date: Date)->Int {
        let timeIntervalInSeconds = date.timeIntervalSince(epoch)
        let forwardDelta = 1/(timeIntervalInSeconds*60*Double(self.fps))
        return Int(ceil(forwardDelta/delta))
    }
    
}
