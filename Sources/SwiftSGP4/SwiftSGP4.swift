import SwiftCelesTrak
import SGP4
import SceneKit

public class SwiftSGP4 {
    private let pi:Double = 3.14159265358979323846
    private let minDelta:Double = 1/60 // currently set as a second
    private let rad:Double
    private var targets:[CelesTrakTarget]
    private var satRecs:[elsetrec]
    private var coordinates:[[SCNVector3]]
    
    public init(_ targets: [CelesTrakTarget], _ date: Date) {
        self.targets = targets
        self.satRecs = [elsetrec]()
        self.rad = 180.0/self.pi
        self.coordinates = [[SCNVector3]]()
    }

    public func propagateOmms(_ targets: [CelesTrakTarget])->[[SCNVector3]] {
        var output = [[SCNVector3]]()
        let count = targets.count
        let today = Date()
        DispatchQueue.concurrentPerform(iterations: count, execute:  { i in
            output.append(computeITRF(targets[i], today, wgs84))
        })
        return output
    }
    
    public func computeITRF(_ target: CelesTrakTarget, _ toDate: Date, _ grabConst:gravconsttype)->[SCNVector3] {
        // orbital set
        var ro = [Double](repeating: 0, count: 3)
        var vo = [Double](repeating: 0, count: 3)
        // time dimension parameters
//        var tsince:Double = 0.0
//        var startmfe:Double = 0.0
//        var stopmfe:Double = 0.0
        
        // struct to pass to sgp4 function
        var satrec = elsetrec()
        // improved algorithm
        let opsMode:CChar = "i".cString(using: .ascii)![0]
        

        // Convert Jd
        let epoch = timestampToJD(target.EPOCH)
//        invjday_SGP4(jd, jdFrac, &year, &month, &day, &hour, &minutes, &seconds)
        
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
        
        // Calculate the target's initial state
        sgp4(&satrec, 0.0, &ro, &vo)
        // transform from TEME to GTRF
        var RGtrf = [Double](repeating: 0, count: 3)
        teme2ecef(&ro, epoch, &RGtrf)
        return [SCNVector3(Float(RGtrf[0]), Float(RGtrf[1]), Float(RGtrf[2]))]
    }
    
    private func timestampToJD( _ dateString: String)->Double {
        let dateFormat = DateFormatter()
        dateFormat.timeZone = TimeZone(abbreviation: "UTC")
        dateFormat.dateFormat = "yyyy-MM-dd'T'HH:mm:ss.SSSSSS"
        let date = dateFormat.date(from: dateString)!

        let calendar = Calendar.current
        let  components  =  calendar.dateComponents([.year, .month, .day, .hour, .second],  from:  date)
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
