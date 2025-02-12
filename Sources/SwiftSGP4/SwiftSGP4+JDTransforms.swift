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

    public func getJulianFromUnix(_ unixSecs: Double) -> Double {
        return (unixSecs / 86400.0) + 2440587.5
    }

    public func getUnixFromJulian(_ julian: Double) -> UInt64 {
        return UInt64((julian - 2440587.5) * 86400.0 + 0.5)
    }

    public func date2jday(_ dt: Date) -> Double {
        let dateUnix = dt.timeIntervalSince1970 * 1000
        let deltaMSec = (dateUnix - unix2000)
        let jd = jd2000 + deltaMSec / 86400000.0
        return jd
    }

    public func jday2date(_ jd: Double) -> Date {
        let jdDeltaMSec = (jd - jd2000) * 86400000.0
        let date = addMilliseconds(utc2000, milliseconds: jdDeltaMSec)
        return date
    }

    public func addDays(_ dt: Date, days: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.day = days
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    public func addMinutes(_ dt: Date, minutes: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.minute = minutes
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    public func addSeconds(_ dt: Date, seconds: Int) -> Date {
        var dateComponents = DateComponents()
        dateComponents.second = seconds
        let updatedDate = calendar.date(byAdding: dateComponents, to: dt)!
        return updatedDate
    }

    public func addMilliseconds(_ dt: Date, milliseconds: Double) -> Date {
        let updatedDate = dt.addingTimeInterval(milliseconds / 1000)
        return updatedDate
    }
    
public func dateString2Date( _ dateString: String)->Date? {
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        return dateFormat.date(from: dateString)
    }
    
    public func timestampToJD( _ date: Date)->Double {
        let calendar = Calendar.current
//        calendar.timeZone = TimeZone(identifier: "UTC")
//        let components = calendar.dateComponents([.year, .month, .day, .hour, .minute, .second],  from:  date)
        let components = calendar.dateComponents(in: TimeZone(identifier: "UTC")!, from: date)
        let year = components.year!.int32()
        let month = components.month!.int32()
        let day = components.day!.int32()
        let hour = components.hour!.int32()
        let minutes = components.minute!.int32()
        let seconds = Double(components.second!)
        let milliseconds = Double(components.nanosecond!)/1000000.0
//        print(date)
//        print("timestamp to JD: \(day)/\(month) \(hour):\(minutes):\(seconds): ms \(milliseconds)")
        var jd:Double = 0.0
        var jdFrac:Double = 0.0
        
        jday(year, month, day, hour, minutes, seconds, &jd, &jdFrac)
        // BUG: add the millisecond to the jdate
        jdFrac = jdFrac + milliseconds/86400000.0
        return jd + jdFrac
    }

extension SwiftSGP4 {
    
    public func framesSinceEpoch( _ date: Date, _ index: Int)->Int {
        print("self.epochs[index]: \(self.epochs[index])")
        let secondsSinceDelta = date.timeIntervalSince(self.epochs[index])
        let forwardDelta = 1/(secondsSinceDelta*Double(self.fps))
        return Int(ceil(forwardDelta/(delta)))
    }
    
    public func steerOffset( _ startDelta: Double, _ index: Int, _ factor: Double = 1)-> (Double, Int) {
        let offsetDelta = startDelta - self.delta*factor*60
        if offsetDelta < 0 {
            return (offsetDelta, index)
        } else {
            return steerOffset(startDelta + delta, index + 1, factor + 1)
        }
    }

    
    
}
