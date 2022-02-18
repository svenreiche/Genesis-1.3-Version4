#ifndef __G4_SEMAFILE_H
#define __G4_SEMAFILE_H

#include <string>

class SemaFile {
public:
	SemaFile();
	~SemaFile();

	void remove(std::string);
	void put(std::string);

private:
	int my_rank_;
};
#endif // __G4_SEMAFILE_H
