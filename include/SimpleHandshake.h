#ifndef G4_SIMPLEHANDSHAKE_H
#define G4_SIMPLEHANDSHAKE_H

class SimpleHandshake {
public:
	SimpleHandshake();
	~SimpleHandshake() = default;
	
	bool doit(const std::string prefix);

private:
	void join(void);
	void drop_file(const std::string fn);
	void wait_for_file(const std::string fname);

	std::string fn_wait_;
	std::string fn_resume_;

	int my_rank_;
	bool busy_wait_ {false};
};

#endif // G4_SIMPLEHANDSHAKE_H
