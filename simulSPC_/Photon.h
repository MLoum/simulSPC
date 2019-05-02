#pragma once
class Photon
{
public:
	__int64 macrotime_;
	unsigned short microtime_;
	unsigned char channel_;
	Photon() {};
	Photon(__int64 macrotime, unsigned short microtime, unsigned char channel) 
	{
		macrotime_ = macrotime;  microtime_ = microtime;  channel_ = channel;
	}
	__int64 compress_data() {return macrotime_ + (microtime_) >> 40 + (channel_ >> 56); }
	~Photon() {};
};

