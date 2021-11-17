#include "PassedEvent.h"
#include "RootIO.h"
#include "Plots.h"
#include "Analysis.h"

#include <vector>
using namespace std;

int main()
{
  //Analysis::RunBigAnalysis();
  vector<PassedEvent> events;
  RootIO::ReadEvents(events);
  //Plots::Plot1(events);
  Plots::Plot2(events);
  //Plot1(events);  
}
