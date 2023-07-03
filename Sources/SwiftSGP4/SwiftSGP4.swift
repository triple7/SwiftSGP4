import SwiftCelesTrak
import SGP4
import simd
import Foundation

public class SwiftSGP4 {
    private let pi:Double = 3.14159265358979323846
    private let rad:Double
    private var targets:[CelesTrakTarget]
    private var satRecs:[elsetrec]
    private lazy var coordinates:[SIMD3<Double>] = {
        return [SIMD3<Double>]()
    }()
    
    // improved algorithm
    private let opsMode:CChar = "i".cString(using: .ascii)![0]
// last time propagation was calculated
    private var lastDate:Date?

    public init(_ targets: [CelesTrakTarget]) {
        self.targets = targets
        self.satRecs = [elsetrec]()
        self.rad = 180.0/self.pi
    }

    public func propagateOmms(_ targets: [CelesTrakTarget], _ secondsFromEpoch: Int, _ fps: Int, _ minDelta: Double = 1/60 /* seconds */)->[[SIMD3<Double>]] {
        let epoch = targets.first!.EPOCH
        lastDate = dateString2Date(epoch)
        // Convert Jd
        let jdEpoch = timestampToJD(epoch)

                                   let count = targets.count
        var output = [[SIMD3<Double>]](repeating: [zeroSimd], count: targets.count)

        // time dimension parameters
        // We are propagating from
        // epoch to secondsFromEpoch by frames per second
        // such that:
        // delta = (1/(60*fps)
        // and count = seconds*fps
        let delta:Double = 1/Double(secondsFromEpoch*fps)
        let dCount = secondsFromEpoch*fps
            
        DispatchQueue.concurrentPerform(iterations: count, execute:  { i in
            output[i] = computeITRF(targets[i], jdEpoch, delta, dCount, wgs84)
        })
        var distances = [Double]()
        for o in output {
            let v = o.first!
            let distance = sqrt(v.x*v.x + v.y*v.y + v.z*v.z)
            distances.append(distance)
        }
        print("Min distance \(distances.min()!)")
        print("max distance: \(distances.max()!)")
        return output
    }
    
    private let zeroSimd = SIMD3<Double>([0, 0, 0])
    public func computeITRF(_ target: CelesTrakTarget, _ epoch: Double, _ delta: Double, _ dCount
                            : Int, _ grabConst:gravconsttype)->[SIMD3<Double>] {
        var output = [SIMD3<Double>](repeating: zeroSimd, count: dCount)

        // struct to pass to sgp4 function
        var satrec = elsetrec()
        
        // populate satrec
//        satrec.classification = target.CLASSIFICATION_TYPE.cString(using: .utf8)![0]
//        let arr = target.OBJECT_NAME.cString(using: .utf8)!
//        satrec.ephtype = target.EPHEMERIS_TYPE.int32()
//        satrec.elnum = target.ELEMENT_SET_NO
//        satrec.revnum = target.REV_AT_EPOCH
        let arr2 = target.OBJECT_ID.cString(using: .utf8)!
        
        var tuple2 = (arr2[0], arr2[1], arr2[2], arr2[3], arr2[4], arr2[5])
//        satrec.satnum = tuple2
        // Initialize sgp4 with the current parameters
        _ = sgp4init(wgs84, opsMode, &tuple2
                 , epoch, target.BSTAR, target.MEAN_MOTION_DOT, target.MEAN_MOTION_DDOT, target.ECCENTRICITY, target.ARG_OF_PERICENTER, target.INCLINATION, target.MEAN_ANOMALY, target.MEAN_MOTION, target.RA_OF_ASC_NODE, &satrec)
        
        // Calculate the target states from epoch to secondsFromEpoch
        DispatchQueue.concurrentPerform(iterations: dCount, execute:  { i in
            // orbital set
            var ro = [Double](repeating: 0, count: 3)
            var vo = [Double](repeating: 0, count: 3)

            let lastSince = Double(i)*delta
//            let epoch1 = epoch + lastSince
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
        
        jday_SGP4(year, month, day, hour, minutes, seconds, &jd, &jdFrac)
        let epoch = jd + jdFrac
return epoch
    }
    
}
