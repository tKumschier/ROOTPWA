#ifndef RANDOMNUMBERGENERATOR_HH_
#define RANDOMNUMBERGENERATOR_HH_

#include<TRandom3.h>

namespace rpwa {

	class randomNumberGenerator {

	  public:

		static randomNumberGenerator* instance();
		TRandom3* getGenerator() { return &_rndGen; }

		unsigned int seed()                     { return _rndGen.GetSeed();     }
		void         setSeed(unsigned int seed) {        _rndGen.SetSeed(seed); }

		double rndm(); // uniform ]0, 1]

	  private:

		randomNumberGenerator() { }
		virtual ~randomNumberGenerator() { }

		static randomNumberGenerator* _randomNumberGenerator;

		TRandom3 _rndGen;

	};

}

#endif
