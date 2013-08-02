#ifndef DATAFORMATS_HCALRECHIT_HCALSOURCEPOSITIONDATA_H
#define DATAFORMATS_HCALRECHIT_HCALSOURCEPOSITIONDATA_H 1

#include <string>

#include "DataFormats/CaloRecHit/interface/CaloRecHit.h"

class HcalSourcePositionData {
public:

  HcalSourcePositionData();
  ~HcalSourcePositionData(){}
  
  inline double indexCounter() const { return indexCounter_; }
  inline double messageCounter() const { return messageCounter_; }
  inline double motorCurrent() const { return motorCurrent_; }
  inline double motorVoltage() const { return motorVoltage_; } 
  inline double reelCounter() const { return reelCounter_; } 
  inline double status() const { return status_; }
  inline double tubeId() const { return tubeId_; }
  inline double driverId() const { return driverId_; }
  inline double sourceId() const { return sourceId_; }
  inline std::string tubeNameFromCoord() const { return tubeNameFromCoord_; }
  inline std::string tubeNameFromSD() const { return tubeNameFromSD_; }
  inline std::string lastCommand() const { return lastCommand_; }
  inline std::string message() const { return message_; }

  void getDriverTimestamp(int& seconds, int& useconds) const;
  void getDAQTimestamp(int& seconds, int& useconds) const;

  void set(double message_counter,
		 int timestamp1_sec,
		 int timestamp1_usec,
		 int timestamp2_sec,
		 int timestamp2_usec,
		 double status,
		 double index_counter,
		 double reel_counter,
		 double motor_current,
		 double motor_voltage,
		 double tube_id,
		 double driver_id,
     double source_id,
     std::string tubeNameFromCoord,
     std::string tubeNameFromSD,
     std::string lastCommand,
     std::string message);

private:
  double messageCounter_;
  double indexCounter_;
  double reelCounter_;
  int timestamp1_sec_;
  int timestamp1_usec_;
  int timestamp2_sec_;
  int timestamp2_usec_;
  double status_;
  double motorCurrent_;
  double motorVoltage_;
  double tubeId_;
  double driverId_;
  double sourceId_;
  std::string tubeNameFromCoord_;
  std::string tubeNameFromSD_;
  std::string lastCommand_;
  std::string message_;
};

std::ostream& operator<<(std::ostream& s, const HcalSourcePositionData& hspd);

#endif
