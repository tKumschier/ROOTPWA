///////////////////////////////////////////////////////////////////////////
//
//    Copyright 2010
//
//    This file is part of rootpwa
//
//    rootpwa is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    rootpwa is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with rootpwa.  If not, see <http://www.gnu.org/licenses/>.
//
///////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------
// File and Version Information:
// $Rev:: 836                         $: revision of last commit
// $Author:: schmeing                 $: author of last commit
// $Date:: 2011-12-21 12:31:38 +0100 #$: date of last commit
//
// Description:
//      Code file for the CompassPwaFilePhaseSpaceIntegrals class that provides
//		functionality to read in phase space integrals from
//		txt files generated by Compass pwa and store them
//
//
// Author List:
//      Stephan Schmeing          TUM            (original author)
//
//
//-------------------------------------------------------------------------

#include <iostream>
#include <sstream>

#include "reportingUtils.hpp"

#include "CompassPwaFilePhaseSpaceIntegrals.h"

using namespace std;
using namespace rpwa;

bool CompassPwaFilePhaseSpaceIntegrals::_Debug = false;

// Default constructor
CompassPwaFilePhaseSpaceIntegrals::CompassPwaFilePhaseSpaceIntegrals(){
}

// Destructor
CompassPwaFilePhaseSpaceIntegrals::~CompassPwaFilePhaseSpaceIntegrals(){
}

// Fills the PhaseSpaceIntegral of the wave with name WaveName into Destination and returns false if WaveName could not be found
bool CompassPwaFilePhaseSpaceIntegrals::PhaseSpaceIntegral( double &Destination, const std::string &WaveName ) const{
	map<const string,const double>::const_iterator it = _PhaseSpaceIntegralsMap.find( WaveName );
	if( it != _PhaseSpaceIntegralsMap.end() ){
		Destination = it->second;
		return true;
	}
	else{
		printErr << "WaveName \"" << WaveName << "\" could not be found in phase space integral for mass bin (" << MassBinStart() << '-' << MassBinEnd() <<")\n";
		return false;
	}
}

// Reads the rest of the information from a phase space integral file stream and returns 0 if no error occurred or a negative number as the error code
bool CompassPwaFilePhaseSpaceIntegrals::ReadIn( std::istream& File ){
	bool Succesful = true; // Is set to false if an error occurs and returned at the end of the function
	stringstream LineStream;
	string Line;

	// Get number of waves
	unsigned int NumWaves = 0;

	if( GetNextValidLine( File, LineStream ) ){
		// Line example between "": "          63"
		LineStream >> NumWaves;

		if( !NumWaves ){
			printErr << "Number of waves either 0 or not an unsigned int\n";
			Succesful = false;
		}
		if( !LineStream.eof() ){
			printWarn << "Number of waves entry longer than expected\n";
		}
	}
	else{
		printErr << "No valid line could be found anymore, but the number of waves was expected\n";
		Succesful = false;
	}

	// Get a map of phase space integrals sorted by their wave names
	char WaveNameCStr[61];
	string WaveName;
	WaveName.reserve(61);
	unsigned int LastNonEmptyCharacter;
	double IntegralValue;
	char apostrophe1; // Takes the first bracket character, which should be a apostrophe
	char apostrophe2; // Takes the second bracket character, which should be a apostrophe

	for( unsigned int i = NumWaves; i--; ){
		if( GetNextValidLine( File, LineStream ) ){
			// Line example between "": "'1-(0-+)0+ rho pi P                                          '  0.890506E-02"
			LineStream >> apostrophe1;
			LineStream.get( WaveNameCStr, 61 );
			LineStream >> apostrophe2 >> IntegralValue;

			if( LineStream.fail() ){
				printErr << "Reading error in a line that is supposed to contain a phase space integral\n";
				Succesful = false;
			}
			else{
				if( !LineStream.eof() ){
					printWarn << "Phase space integral line longer than expected\n";
				}
				if( apostrophe1 != '\'' || apostrophe2 != '\'' ){
					printWarn << "Wave name bracket not an apostrophe\n";
					if( _Debug ){
						printDebug << "Brackets: '" << apostrophe1 << "','"<< apostrophe2 << "'\n";
					}
				}

				WaveName = WaveNameCStr;

				// Chop of following whitespaces
				LastNonEmptyCharacter = WaveName.find_last_not_of(' ');
				WaveName.resize(LastNonEmptyCharacter + 1);

				if( !( _PhaseSpaceIntegralsMap.insert( pair<const string, const double>(WaveName, IntegralValue) ).second ) ){
					printErr << "Two waves with the same name in the file: " << WaveName << '\n';
					Succesful = false;
				}
			}
		}
		else{
			printErr << "Less wave names found than specified\n";
			Succesful = false;
		}
	}

	// If there are still valid lines, the number of waves is not correct
	if( GetNextValidLine( File, Line ) ){
		// Some editors put automatically an empty line at the end, so this error has to be caught
		if( !Line.empty() || GetNextValidLine( File, Line ) ){
			printErr << "More wave names found than specified \n";
			Succesful = false;
		}
	}

	return Succesful; // No error occurred
}

// Prints all important variables of class
ostream& CompassPwaFilePhaseSpaceIntegrals::Print( ostream& Out ) const{
	CompassPwaFileBase::Print( Out );

	Out << "Number of waves: " << _PhaseSpaceIntegralsMap.size() << '\n';

	for ( map<const string,const double>::const_iterator it = _PhaseSpaceIntegralsMap.begin() ; it != _PhaseSpaceIntegralsMap.end(); it++ ){
		Out << '\'' << it->first << '\'' << ' ' << it->second << '\n';
	}

	return Out;
}

// Combines the matching integrals from Integrals to one for the given mass bin and given waves, stores it in Destination and returns a reference to Destination
bool CompassPwaFilePhaseSpaceIntegrals::Combine( vector<double>& Destination, const deque<const CompassPwaFilePhaseSpaceIntegrals *>& Integrals, const vector< vector<string> >& WaveNames, double MassBinStart, double MassBinEnd ){
	bool Succesful = true;

	// Determines total number of waves
	unsigned int TotNumWaves = 0;
	for( unsigned int s=0; s < WaveNames.size(); ++s ){
		TotNumWaves += WaveNames[s].size();
	}

	Destination.clear();
	Destination.resize( TotNumWaves );

	// For now it just returns the integral in the middle of the mass bin specified
	unsigned int middle = Integrals.size() / 2;

	unsigned int j=0;
	for( unsigned int s = 0; s < WaveNames.size(); ++s ){ // Section
		for( unsigned int i = 0; i < WaveNames[s].size(); ++i ){
			Succesful = Succesful && Integrals[middle]->PhaseSpaceIntegral( Destination[j++], WaveNames[s][i] );
		}
	}

	return Succesful;
}