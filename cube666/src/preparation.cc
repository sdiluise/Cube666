#include "preparation.hh"
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "G4Poisson.hh"
#include "istream"




G4ThreeVector randomPos2D(G4double radius,G4double z){
  G4double phi = 2*TMath::Pi()*drand48();
  G4double r   = radius*sqrt(drand48());
  G4ThreeVector vec(r*cos(phi),r*sin(phi),z);
  return vec;
}




G4ThreeVector sphericalDistribution(const char* direction){
  G4double cosTheta = 2*drand48() - 1;
  G4double sinTheta = sqrt(1-cosTheta*cosTheta);
  G4double phi = 2*PI*drand48();

  G4ThreeVector vector;

  if(!strcmp(direction,"negative"))
    vector.set(sinTheta*cos(phi),sinTheta*sin(phi),-abs(cosTheta));
  else if(!strcmp(direction,"positive"))
    vector.set(sinTheta*cos(phi),sinTheta*sin(phi),abs(cosTheta));
  else vector.set(sinTheta*cos(phi),sinTheta*sin(phi),cosTheta);

  return vector;
}




G4ThreeVector polarizationVec(G4ThreeVector momDir){
  //generate a random 3D vector perpendicular to momDir

  if(!momDir.mag2()) return sphericalDistribution();

  G4ThreeVector polarization = momDir.orthogonal();
  G4ThreeVector perp         = momDir.cross(polarization);

  G4double phi = twopi*G4UniformRand();
  polarization = std::cos(phi)*polarization + std::sin(phi)*perp;

  return polarization.unit();
}







vector<G4ThreeVector> setGridWireVector(const char* xyAxis,vector<double>& wireHalfLength,
					double plateInnerR, double wirePitch, double wireOuterR){


  //int nwires = (int)((CATHODE_PLATE_INNER_RADIUS-CATHODE_WIRE_PITCH/2)/CATHODE_WIRE_PITCH)+1;
  //double R = CATHODE_PLATE_INNER_RADIUS;


  int nwires = (int)((plateInnerR-wirePitch/2)/wirePitch)+1;
  double R = plateInnerR;

  vector<G4ThreeVector> pos;

  wireHalfLength.clear();

  double z;
  
  if(!strcmp(xyAxis,"x"))      z = 0;
  else if(!strcmp(xyAxis,"y")) z = 2*wireOuterR+.1*mm;
  for(int i=-nwires;i<nwires;i++){
    double coord = (!strcmp(xyAxis,"y") ? CUBE666_CATHODE_WIRE_PITCH_Y : CUBE666_CATHODE_WIRE_PITCH_X)  *(i+1./2);
    double halflength = sqrt(R*R - coord*coord);

    if(!halflength) continue;

    if(!strcmp(xyAxis,"x"))      pos.push_back(G4ThreeVector(coord,0.,z));
    else if(!strcmp(xyAxis,"y")) pos.push_back(G4ThreeVector(0.,coord,z));

    wireHalfLength.push_back(halflength);
  }


  return pos;
}










double getWLSThickness(double convEff){

  return ((convEff) > 0 && (convEff) < 1)?(-TMath::Log(1-(convEff))*(WLS_MEAN_ABSORPTION_LENGTH)):(WLS_THICKNESS_100_PERCENT_CONVERSION_EFFICIENCY);

}







vector<TVector2> setPMTVector2D(){

  int nx = CUBE666_NPMT_X;
  int ny = CUBE666_NPMT_Y;

  double pitchx = CUBE666_PMT_PITCH_X;
  double pitchy = CUBE666_PMT_PITCH_Y;

  double offsetx = nx%2+1;
  double offsety = ny%2+1;


  vector<TVector2> pos;

  for(int xi=-(nx/2-offsetx);xi<=nx/2;xi++){
    for(int yi=-(ny/2-offsety);yi<=ny/2;yi++){
      double posx = (xi-offsetx/2)*pitchx;
      double posy = (yi-offsety/2)*pitchy;
      pos.push_back(TVector2(posx,posy));
    }
  }

  return pos;
}





vector<G4ThreeVector> setPMTVector(){

  vector<G4ThreeVector> rPMT;

  G4double posz = CUBE666_PMT_POS_Z_LARCOL;
  
  vector<TVector2> rPMT2D = setPMTVector2D();

  for(int i=0;i<rPMT2D.size();i++){
    rPMT.push_back(G4ThreeVector(rPMT2D[i].X(),rPMT2D[i].Y(),posz));
  }

  return rPMT;
}






/*
vector<vector<string> > readTextfile_string(string filename,char delimiter){
  //hrm, somehow this function doesn't work.
  //--> go figure !!


  ifstream inputFile(filename.c_str());
  
  vector<vector<string> > output;
  
  if(!inputFile) return output;

  string line;
  while(getline(inputFile,line)){
    if(line.empty() || line == '\r\n' || line == '\n'  || line == '\r' || line.find("#") != string::npos) continue;
    //if(line.empty() || line == '\r\n' || line == '\n'  || line == '\r') continue;
    vector<string> row;
    string value;
    istringstream iss(line,istringstream::in);

    if(delimiter == '\0'){
      while(!iss.eof()){
	iss >> value;
	row.push_back(value);
      }

    }else{
      while(getline(iss,value,delimiter)) row.push_back(value);
    }

    output.push_back(row);
  }


  return output;
}


*/


/*
vector<vector<double> > readTextfile_float(string filename,char delimiter){

  ifstream inputFile(filename.c_str());
  
  vector<vector<double> > output;
  
  if(!inputFile) return output;
  
  string line;
  int count=0;
  while(getline(inputFile,line)){
    if(line.empty() || line == '\r\n' || line == '\n'  || line == '\r' || line.find("#") != string::npos) continue;
    vector<double> row;
    string value;
    istringstream iss(line,istringstream::in);

    if(delimiter == '\0'){
      while(!iss.eof()){
	iss >> value;
	row.push_back(atof(value.c_str()));
      }

    }else{
      while(getline(iss,value,delimiter)) row.push_back(atof(value.c_str()) );
    }
    
    output.push_back(row);

  }
  
  return output;
}

*/
