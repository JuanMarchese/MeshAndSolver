/*
 * ConfigurationReader.cpp
 *
 *  Created on: Nov 2, 2014
 *      Author: juan
 */

#include "ConfigurationReader.h"

ConfigurationReader::ConfigurationReader() {
	this->commentCharacter = '#';
	this->separationCharacter = '=';
}

ConfigurationReader::~ConfigurationReader() {
}

void ConfigurationReader::init(const std::string& path,const char& commentCharacter,const char& separationCharacter){
	this->path = path;
	this->commentCharacter = commentCharacter;
	this->separationCharacter = separationCharacter;
}

void ConfigurationReader::process(){

	std::ifstream myfile(this->path);
	std::string line;

	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			processLine(line);
		}
		myfile.close();
	}

}

std::string ConfigurationReader::getParameter(const std::string& name, const std::string& defaultValue){
	std::map<std::string,std::string>::iterator it;
	std::string retVal = defaultValue;
	it = this->parameters.find(name);

	if(it != this->parameters.end()){
		retVal = it->second;
	}

	return retVal;
}

void ConfigurationReader::processLine(const std::string& line){
	std::string cleanedLine;
	size_t commentPos = line.find(this->commentCharacter);

	cleanedLine = line.substr(0,commentPos);

	size_t spacePos = cleanedLine.find(' ');

	while(spacePos != cleanedLine.npos){
		cleanedLine.erase(spacePos,1);
		spacePos = cleanedLine.find(' ');
	}

	size_t separatorPos = cleanedLine.find(this->separationCharacter);
	if(separatorPos != cleanedLine.npos){
		std::string key = cleanedLine.substr(0,separatorPos);
		std::string value = cleanedLine.substr(separatorPos+1,cleanedLine.npos);
		parameters[key] = value;
	}

}
