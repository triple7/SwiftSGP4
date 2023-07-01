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

    public func propagateOmms(_ targets: [CelesTrakTarget], _ secondsFromEpoch: Int, _ fps: Int, _ minDelta: Double = 1/60 /* seconds */)->[SIMD3<Double>] {
        var output = [SIMD3<Double>]()
        let count = targets.count
        DispatchQueue.concurrentPerform(iterations: count, execute:  { i in
            output.append(contentsOf: computeITRF(targets[i], secondsFromEpoch, fps, wgs84))
        })
        return output
    }
    
    public func computeITRF(_ target: CelesTrakTarget, _ secondsFromEpoch: Int, _ fps: Int, _ grabConst:gravconsttype)->[SIMD3<Double>] {
        var output = [SIMD3<Double>]()
        // time dimension parameters
        // We are propagating from
        // epoch to secondsFromEpoch by frames per second
        // such that:
        // delta = (1/(60*fps)
        // and count = seconds*fps
        let delta:Double = 1/Double(secondsFromEpoch*fps)
        let dCount = secondsFromEpoch*fps
        // struct to pass to sgp4 function
        var satrec = elsetrec()
        

        // Convert Jd
        let epoch = timestampToJD(target.EPOCH)
        
        // populate satrec
        satrec.classification = target.CLASSIFICATION_TYPE.cString(using: .utf8)![0]
        let arr = target.OBJECT_ID.cString(using: .ascii)!
        let tuple = (arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7], arr[8], arr[9], arr[10])
        satrec.intldesg = tuple
        satrec.ephtype = target.EPHEMERIS_TYPE.int32()
        satrec.elnum = target.ELEMENT_SET_NO
        satrec.revnum = target.REV_AT_EPOCH
        let arr2 = target.OBJECT_ID.cString(using: .utf8)!
        var tuple2 = (arr2[0], arr2[1], arr2[2], arr2[3], arr2[4], arr2[5])
        satrec.satnum = tuple2
        // Initialize sgp4 with the current parameters
        _ = sgp4init(wgs84, opsMode, &tuple2
                 , epoch, target.BSTAR, target.MEAN_MOTION_DOT, target.MEAN_MOTION_DDOT, target.ECCENTRICITY, target.ARG_OF_PERICENTER, target.INCLINATION, target.MEAN_ANOMALY, target.MEAN_MOTION, target.RA_OF_ASC_NODE, &satrec)
        
        // Calculate the target states from epoch to secondsFromEpoch
        DispatchQueue.concurrentPerform(iterations: dCount, execute:  { i in
            // orbital set
            var ro = [Double](repeating: 0, count: 3)
            var vo = [Double](repeating: 0, count: 3)

            let deltaFromEpoch = Double(i)*delta
            sgp4(&satrec, deltaFromEpoch, &ro, &vo)
            // transform from TEME to GTRF
            var RGtrf = [Double](repeating: 0, count: 3)
            teme2ecef(&ro, epoch+deltaFromEpoch, &RGtrf)
            output.append(SIMD3<Double>(RGtrf))
        })
        return output
    }
    
    private func timestampToJD( _ dateString: String)->Double {
        print("Date is \(dateString)")
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        lastDate = dateFormat.date(from: dateString)!

        let calendar = Calendar.current
        let  components  =  calendar.dateComponents([.year, .month, .day, .hour, .second],  from:  lastDate!)
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
