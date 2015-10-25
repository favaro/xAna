#include "interface/setupReader.hh"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

using namespace std;

std::string itostr(int i ) {
  std::stringstream out;
  out << i;
  return out.str();
}

SetupReader::SetupReader( string setupfile ) {
  _setupfile = setupfile;
}

int   SetupReader::getInt(  string vname, bool &found ) const {
  ifstream input(_setupfile.c_str());
  
  while ( input.good() && !input.eof() ) {
    string line;
    getline(input,line,'\n');
    if( !input.good() ) cerr << " SetupReader[ERROR]: input file " << _setupfile << " is not opened correctly (or may not exist)..." << endl;
    /// comment
    if( line.find("#") != string::npos ) continue;
    int posSeparator = line.find("=");    
    if( posSeparator < 0 ) posSeparator = line.find(":");    
    if( posSeparator < 0 ) continue;
    
    
    /// split the line in key : val
    string key = line.substr(0,posSeparator);
    string val = line.substr(posSeparator+1,line.size());
    removeSpaces(key);
    removeSpaces(val);
    
    if( key.find( vname ) != string::npos ) {
      found = true;
      input.close();
      return atoi(val.c_str());
    }
  }
  input.close();
  found = false;
  return -1;
}

float SetupReader::getFloat(  string vname, bool &found ) const {
  ifstream input(_setupfile.c_str());
  
  while ( input.good() && !input.eof() ) {
    string line;
    getline(input,line,'\n');
    if( !input.good() ) cerr << " SetupReader[ERROR]: input file " << _setupfile << " is not opened correctly (or may not exist)..." << endl;

    /// comment
    if( line.find("#") != string::npos ) continue;
    int posSeparator = line.find("=");    
    if( posSeparator < 0 ) posSeparator = line.find(":");    
    if( posSeparator < 0 ) continue;
    
    
    /// split the line in key : val
    string key = line.substr(0,posSeparator);
    string val = line.substr(posSeparator+1,line.size());
    removeSpaces(key);
    removeSpaces(val);
    
    if( key.find( vname ) != string::npos ) {
      found = true;
      input.close();
      return atof(val.c_str());
    }
  }
  input.close();
  found = false;
  return -1.0;
}


string   SetupReader::getString(  string vname, bool &found ) const {
  ifstream input(_setupfile.c_str());
  if( !input.good() ) cerr << " SetupReader[ERROR]: input file " << _setupfile << " is not opened correctly (or may not exist)..." << endl;
  while ( input.good() && !input.eof() ) {
    string line;
    getline(input,line,'\n');
    
    /// comment
    if( line.find("#") != string::npos ) continue;
    int posSeparator = line.find("=");    
    if( posSeparator < 0 ) posSeparator = line.find(":");    
    if( posSeparator < 0 ) continue;
    
    
    /// split the line in key : val
    string key = line.substr(0,posSeparator);
    string val = line.substr(posSeparator+1,line.size());
    removeSpaces(key);
    removeSpaces(val);
    
    if( key.find( vname ) != string::npos ) {
      found = true;
      input.close();
      return val;
    }
  }
  input.close();
  found = false;
  return "default";
}


SetupReader::~SetupReader(void) {}


void SetupReader::removeSpaces(string &str) const {
  unsigned long ipos = -1;
  while( (ipos = str.find(" ")) != string::npos ) str.replace(ipos,1,"");
}

