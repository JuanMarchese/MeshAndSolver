/*
 * ConfigurationReader.h
 *
 *  Created on: Nov 2, 2014
 *      Author: juan
 */

#ifndef CONFIGURATIONREADER_H_
#define CONFIGURATIONREADER_H_

#include <map>
#include <string>
#include <fstream>

class ConfigurationReader {

private:
	std::map<std::string,std::string> parameters;
	std::string path;
	char commentCharacter;
	char separationCharacter;


public:
	ConfigurationReader();
	virtual ~ConfigurationReader();

	void init(const std::string& path,const char& commentCharacter,const char& separationCharacter);
	void process();
	std::string getParameter(const std::string& name, const std::string& defaultValue = "");

private:

	void processLine(const std::string& line);

};

#endif /* CONFIGURATIONREADER_H_ */
