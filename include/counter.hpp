#ifndef COUNTER_HPP
#define COUNTER_HPP



class counter
{
private :
	int step_max;
	int step;
public :
	void set_range(const int &N) {step_max = N;}
	counter(int s_max);
	void reset(void) {step = 0;}
	int show(void) {return step;}
	bool up(void);
	bool good(void) {return (step<step_max) ? true : false;}
	void operator ++() {step++;}

};



#endif // ende der Counter Classe
