#ifndef G4_SIMPLEHANDSHAKE_H
#define G4_SIMPLEHANDSHAKE_H

class SimpleHandshake {
public:
	SimpleHandshake();
	~SimpleHandshake() = default;
	
	bool doit(const std::string prefix);

private:
	void drop_file(const std::string fn);
	void wait_for_file(const std::string fname);

	int my_rank_;
	std::string fn_wait_;
	std::string fn_resume_;
};

#endif // G4_SIMPLEHANDSHAKE_H
