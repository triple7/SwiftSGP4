import SwiftCelesTrak
import SGP4
import simd
import Foundation

public class SwiftSGP4 {
    private let pi:Double = 3.14159265358979323846
    private let rad:Double
    private let deg2rad:Double
    private let xpdotp:Double
    private let xpdotInv:Double
    private let xpdotInv2:Double
    private let minPDay:Double
    private let secPDay:Double

    public var targets:[CelesTrakTarget]
    private var satRecs:[elsetrec]
    public var coordinates:[[SIMD3<Double>]]
    public var lastT:Double = 0
    // improved algorithm
    private let opsMode:CChar = "i".cString(using: .utf8)![0]
// last minute since value time propagation was calculated
    private var lastTSince:Double?

    public init(_ targets: [CelesTrakTarget]) {
        self.targets = targets
        self.satRecs = [elsetrec](repeating: elsetrec(), count: targets.count)
        self.rad = 180.0/self.pi
        self.deg2rad = pi / 180.0
        self.xpdotp = 1440.0/(2.0*pi)
        self.xpdotInv = self.xpdotp*1440
        self.xpdotInv2 = self.xpdotp*1440*1440
        self.minPDay = 1440
        self.secPDay = 1440*60
        self.coordinates = [[SIMD3<Double>]]()
    }

    public func propagateOmms( _ secondsFromEpoch: Int, _ fps: Int, _ minDelta: Double = 1/60 /* seconds */) {
        let epoch = targets.first!.EPOCH
        // Convert Jd
        let jdEpoch = timestampToJD(epoch)

                                   let count = targets.count
        coordinates = [[SIMD3<Double>]](repeating: [zeroSimd], count: targets.count)
        // time dimension parameters
        // We are propagating from
        // epoch to secondsFromEpoch by frames per second
        // such that:
        // delta = (1/(60*fps)
        // and count = seconds*fps
        let delta:Double = 1/Double(secondsFromEpoch*fps)
        let dCount = secondsFromEpoch*fps
            
        DispatchQueue.concurrentPerform(iterations: count, execute:  { i in
            coordinates[i] = computeITRF(i, targets[i], jdEpoch, delta, dCount, wgs84)
        })
    }
    
    private let zeroSimd = SIMD3<Double>([0, 0, 0])
    // Using generic number as we don't need satNum for propagation
    private var genSatNum = (Int8(77), Int8(77), Int8(77), Int8(77), Int8(77), Int8(77))
    public func computeITRF(_ satrecIndex: Int, _ target: CelesTrakTarget, _ epoch: Double, _ delta: Double, _ dCount
                            : Int, _ grabConst:gravconsttype)->[SIMD3<Double>] {
        var output = [SIMD3<Double>](repeating: zeroSimd, count: dCount)

        // struct to pass to sgp4 function
        var satrec = satRecs[satrecIndex]
            // save the last frame from epoch for the next cycle
            lastT = Double(dCount)*delta
            
        // no need to populate satrec
        // But this is how properties are transformed
//        satrec.classification = target.CLASSIFICATION_TYPE.cString(using: .utf8)![0]
//        let arr = target.OBJECT_NAME.cString(using: .utf8)!
//        satrec.ephtype = target.EPHEMERIS_TYPE.int32()
//        satrec.elnum = target.ELEMENT_SET_NO
//        satrec.revnum = target.REV_AT_EPOCH
//        let arr2 = target.OBJECT_ID.cString(using: .utf8)!
        
//        var tuple2 = (arr2[0], arr2[1], arr2[2], arr2[3], arr2[4], arr2[5])
        // Initialize sgp4 with the current parameters
        _ = sgp4init(wgs72, opsMode, &genSatNum
                 , epoch, target.BSTAR, target.MEAN_MOTION_DOT/xpdotInv, target.MEAN_MOTION_DDOT/xpdotInv2, target.ECCENTRICITY*deg2rad, target.ARG_OF_PERICENTER*deg2rad, target.INCLINATION*deg2rad, target.MEAN_ANOMALY*deg2rad, target.MEAN_MOTION/xpdotp, target.RA_OF_ASC_NODE*deg2rad, &satrec)
        
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
            output[i] = SIMD3<Double>(RGtrf)
        })
        return output
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
