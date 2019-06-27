
#ifndef RSMDATA_H_
#define RSMDATA_H_


#define MAXL  (2000)
#define MAXL0 (2001)

class rsmData{

public:
  rsmData();
  ~rsmData();

private:
  char status[20];
  float bprob;
  int count;
  int m0status;
  int m1status;
  char mActive[MAXL0];
  char mPassive[MAXL0];
  int errcode;
  bool deterministicBind;
  bool deterministicExec;
  char product[MAXL0];

};


#endif
