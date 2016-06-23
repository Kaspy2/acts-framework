//
//  DataClassOne
//  ACTFW
//
//  Created by Andreas Salzburger on 11/05/16.
//
//

#ifndef DataClassOne_h
#define DataClassOne_h

#include <memory>
#include <string>
#include <sstream>
#include <vector>

namespace FWE {
    
    /// @class Simple data class with 
    class DataClassOne {
        
    public :
        /// Constructor
        DataClassOne(const std::string& stringData, size_t eventData) :
          m_dataString(stringData),
          m_dataSizeT(eventData)
        {}
        
        /// Destructor
        ~DataClassOne(){}
        
        /// the contained data : string 
        const std::string data() const;

        
    private:
        std::string m_dataString; ///< data member string
        size_t      m_dataSizeT;  ///< data member size_t
        
    };
    
    inline const std::string DataClassOne::data() const
    {
        std::ostringstream oss;
        oss << "Data : " << m_dataString << " | " << m_dataSizeT;
        return oss.str();
    }
    
    typedef std::vector< std::unique_ptr<DataClassOne> > DataClassOneCollection;
    
}


#endif