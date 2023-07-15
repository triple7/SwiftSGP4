import SwiftCelesTrak
import SGP4
import simd
import Foundation

public class SwiftSGP4 {
    private let pi:Double = 3.14159265358979323846
    private let jd1950:Double = 2433281.5
    private let rad:Double
    private let deg2rad:Double
    private let xpdotp:Double
    private let xpdotInv:Double
    private let xpdotInv2:Double
    private let minPDay:Double
    private let secPDay:Double

    public var targets:[CelesTrakTarget]
    private var satRecs:[elsetrec]
    public var coordinates:ContiguousArray<ContiguousArray<SIMD3<Double>>>
    // improved algorithm
    private let opsMode:CChar = "a".cString(using: .utf8)![0]
// last minute since value time propagation was calculated
    private var lastTSince:Double = 0
    // default second since last and fps
    private let secondsFromEpoch:Int = 1
    private let fps:Int = 30

    private let targetCount:Int
    private let epoch:String
    private let jdEpoch:Double

    
    public init(_ targets: [CelesTrakTarget]) {
        
        self.targets = targets
        self.targetCount = targets.count
        self.rad = 180.0/self.pi
        self.deg2rad = pi / 180.0
        self.minPDay = 1440
        self.secPDay = 1440*60
        self.xpdotp = self.minPDay/(2.0*pi)
        self.xpdotInv = self.xpdotp*self.minPDay
        self.xpdotInv2 = self.xpdotp*self.minPDay*self.minPDay

 epoch = targets.first!.EPOCH
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        let date = dateFormat.date(from: epoch)!

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
        jdEpoch = jd + jdFrac
        print("JD epoch: \(jdEpoch)")
        print("Btar: \(targets[0].BSTAR)")

        self.satRecs = [elsetrec]()
        var target = targets.first!
        print("jd \(jd)")
        print(" jdfrac: \(jdFrac)")
        print("bstar \(target.BSTAR)")
        print("mean motion dot\(target.MEAN_MOTION_DOT)")
        print("mean motion ddot\(target.MEAN_MOTION_DDOT)")
        print("eccentricity\(target.ECCENTRICITY)")
        print("argument of pericenter\(target.ARG_OF_PERICENTER)")
        print("Inclination \(target.INCLINATION)")
        print("mean anomaly \(target.MEAN_ANOMALY)")
        print("mean motion\(target.MEAN_MOTION)")
        print("rra of node\(target.RA_OF_ASC_NODE)")

        for target in targets {
            var satrec = elsetrec()
            satrec.elnum = target.ELEMENT_SET_NO
            satrec.revnum = target.REV_AT_EPOCH
            satrec.classification = target.CLASSIFICATION_TYPE.cString(using: .unicode)![0]
            
            _ = sgp4init(wgs72, opsMode, &genSatNum
                         , jdEpoch - jd1950, target.BSTAR, target.MEAN_MOTION_DOT/xpdotInv, target.MEAN_MOTION_DDOT/xpdotInv2, target.ECCENTRICITY, target.ARG_OF_PERICENTER*deg2rad, target.INCLINATION*deg2rad, target.MEAN_ANOMALY*deg2rad,
                         target.MEAN_MOTION/xpdotp, target.RA_OF_ASC_NODE*deg2rad, &satrec)

            
            satRecs.append(satrec)
        }

        let targetFrames = ContiguousArray<SIMD3<Double>>(repeating: zeroSimd, count: secondsFromEpoch*fps)
        self.coordinates = ContiguousArray<ContiguousArray<SIMD3<Double>>>(repeating: targetFrames, count: targetCount)
    }

    public func propagateOmms( _ minDelta: Double = 1/60 /* seconds */) {
        
                                   let count = targets.count
        // time dimension parameters
        // We are propagating from
        // epoch to secondsFromEpoch by frames per second
        // such that:
        // delta = (1/(60*fps)
        // and count = seconds*fps
        let delta:Double = 1/Double(secondsFromEpoch*fps)
        let dCount = secondsFromEpoch*fps
        // save the last frame from epoch for the next cycle
        lastTSince = Double(dCount)*delta

            
        DispatchQueue.concurrentPerform(iterations: count, execute:  { i in
            computeITRF(i, jdEpoch, delta, dCount)
        })
    }
    
    private let zeroSimd = SIMD3<Double>([0, 0, 0])
    // Using generic number as we don't need satNum for propagation
    private var genSatNum = (Int8(77), Int8(77), Int8(77), Int8(77), Int8(77), Int8(77))
    public func computeITRF(_ satrecIndex: Int, _ epoch: Double, _ delta: Double, _ dCount
                            : Int) {

        // struct to pass to sgp4 function
        var satrec = satRecs[satrecIndex]
            
        // Calculate the target states from epoch to secondsFromEpoch
        DispatchQueue.concurrentPerform(iterations: dCount, execute:  { i in
            // orbital set
            var ro = [Double](repeating: 0, count: 3)
            var vo = [Double](repeating: 0, count: 3)

            let lastSince = Double(i)*delta
            sgp4(&satrec, lastSince, &ro, &vo)
            // transform from TEME to GTRF
            var RGtrf = [Double](repeating: 0, count: 3)
            let gmst = gstime(jdut1: epoch)
            let gmstCos = cos(gmst)
            let gmstSin = sin(gmst)

            teme2ecefOptimised(&ro, epoch, gmstCos, gmstSin, &RGtrf)
            self.coordinates[satrecIndex][i] = SIMD3<Double>(RGtrf)
        })
//        if targets[satrecIndex].NORAD_CAT_ID == 25544 {
//            print("ISS z: \(self.coordinates[satrecIndex][0])")
//        }
    }

    private func dateString2Date( _ dateString: String)->Date {
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        let date = dateFormat.date(from: dateString)!
        return date
    }
    
    private func timestampToJD( _ dateString: String)->Double {
        let date = dateString2Date( dateString)
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
        let epoch = jd + jdFrac
return epoch
    }
    
}
