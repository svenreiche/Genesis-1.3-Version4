#ifndef G4_SIMPLEHANDSHAKE_H
#define G4_SIMPLEHANDSHAKE_H

#include <map>
#include <string>

class SimpleHandshake {
public:
	SimpleHandshake();
	~SimpleHandshake() = default;
	
	bool doit(std::map<std::string,std::string> *, const std::string prefix);

private:
	void usage(void);
	void join(void);
	void drop_file(const std::string fn);
	void wait_for_file(const std::string fname);

	std::string ext_wait_;
	std::string ext_resume_;

	int my_rank_;
	bool busy_wait_ {false};
};

#endif // G4_SIMPLEHANDSHAKE_H
