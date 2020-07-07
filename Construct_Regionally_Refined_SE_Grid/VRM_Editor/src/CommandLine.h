///////////////////////////////////////////////////////////////////////////////
///
///	\file    CommandLine.h
///	\author  Paul Ullrich
///	\version February 24, 2013
///
///	<remarks>
///		Copyright 2000-2010 Paul Ullrich
///
///		This file is distributed as part of the Tempest source code package.
///		Permission is granted to use, copy, modify and distribute this
///		source code and its documentation under the terms of the GNU General
///		Public License.  This software is provided "as is" without express
///		or implied warranty.
///	</remarks>

#ifndef _COMMANDLINE_H_
#define _COMMANDLINE_H_

#include "SQuadGen/Exception.h"

#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>

//=====================================================================================

//  Types of parameters.
//========================
enum ParameterType 
{
    ParameterTypeNone,
    ParameterTypeBool,
    ParameterTypeString,
    ParameterTypeInt,
    ParameterTypeDouble,
};


//=====================================================================================
class CommandLineParameter 
  //  A command line parameter.
  //===============================
{
    public: 
      //  Default constructor.
      //======================
      CommandLineParameter(std::string strName, 
                           std::string strDescription) : m_strName(std::string("--") + strName),
		                                         m_strDescription(strDescription)
      { }

      //  Virtual destructor.
      //=====================
      virtual ~CommandLineParameter() { }

      //  Identify the type of parameter.
      //=================================
      virtual ParameterType GetParameterType() 
      {
        return ParameterTypeNone;
      }

      //  Number of values required.
      //=============================
      virtual int GetValueCount() const = 0;

      //  Print the usage information of this parameter.
      //================================================
      virtual void PrintUsage() const = 0;

      //  Activate this parameter.
      //===========================
      virtual void Activate() { }

      //  Set the value from a string.
      //==============================
      virtual void SetValue(int ix, std::string strValue) 
      {
        _EXCEPTIONT("Invalid value index.");
      }

    public:
      std::string m_strName;          //  Name of this parameter.
      std::string m_strDescription;   //  Description of this parameter.
};
//=====================================================================================


//=====================================================================================
class CommandLineParameterBool : public CommandLineParameter 
  //  A command line boolean.
  //===============================
{
    public:
      //  Constructor.
      //===============
      CommandLineParameterBool(bool & ref, 
                               std::string strName, 
                               std::string strDescription) : CommandLineParameter(strName,strDescription), 
                                                              m_fValue(ref)
      {
        m_fValue = false;
      }

      //  Identify the type of parameter.
      //=================================
      virtual ParameterType GetParameterType() { return ParameterTypeBool; }

      //  Number of values required.
      //============================
      virtual int GetValueCount() const { return (0); }

      //  Print the usage information of this parameter.
      //=================================================
      virtual void PrintUsage() const 
      {
        printf("  %s <bool> [%s] %s\n", m_strName.c_str(), (m_fValue)?("true"):("false"), 
                                        m_strDescription.c_str());
      }

      //  Activate this parameter.
      //===========================
      virtual void Activate() { m_fValue = true; }

    public:
      bool & m_fValue;   //  Parameter value.
};
//=====================================================================================


//=====================================================================================
class CommandLineParameterString : public CommandLineParameter 
  //  A command line string.
  //========================
{
    public:
      //  Constructor.
      //===============
      CommandLineParameterString(std::string & ref,
                                 std::string strName,
                                 std::string strDefaultValue,
                                 std::string strDescription) : CommandLineParameter(strName, strDescription),
                                                               m_strValue(ref)
      {
        m_strValue = strDefaultValue;
      }

      //  Identify the type of parameter.
      //===================================
      virtual ParameterType GetParameterType() { return ParameterTypeString; }

      //  Number of values required.
      //==============================
      virtual int GetValueCount() const { return (1); }

      //  Print the usage information of this parameter.
      //==================================================
      virtual void PrintUsage() const 
      {
        printf("  %s <string> [\"%s\"] %s\n", m_strName.c_str(), m_strValue.c_str(), 
                                              m_strDescription.c_str());
      }

      //  Set the value from a string.
      //===============================
      virtual void SetValue( int ix, std::string strValue) 
      {
        if (ix != 0) { _EXCEPTIONT("Invalid value index."); }
        m_strValue = strValue.c_str();
      }

    public:
      std::string & m_strValue;   //  Parameter value.
};
//=====================================================================================


//=====================================================================================
class CommandLineParameterInt : public CommandLineParameter 
  //  A command line integer.
  //===========================
{
    public:
      //  Constructor.
      CommandLineParameterInt(int &       ref,
                              std::string strName,
                              int         dDefaultValue,
                              std::string strDescription) : CommandLineParameter(strName, strDescription),
                                                             m_dValue(ref)
      {
        m_dValue = dDefaultValue;
      }

      //  Identify the type of parameter.
      //==================================
      virtual ParameterType GetParameterType() const { return ParameterTypeInt; }

      //  Number of values required.
      //==============================
      virtual int GetValueCount() const { return (1); }

      //  Print the usage information of this parameter.
      //==================================================
      virtual void PrintUsage() const 
      {
        printf("  %s <integer> [%i] %s\n", m_strName.c_str(), m_dValue, 
                                           m_strDescription.c_str());
      }

      //  Set the value from a string.
      //================================
      virtual void SetValue( int ix, std::string strValue) 
      {
        if (ix != 0) { _EXCEPTIONT("Invalid value index."); }
        m_dValue = atoi(strValue.c_str());
      }

    public:
      int & m_dValue;   //  Parameter value.
};
//=====================================================================================


//=====================================================================================
class CommandLineParameterDouble : public CommandLineParameter 
  //  A command line double.
  //==========================
{
    public:
      //  Constructor.
      //===============
      CommandLineParameterDouble(double &    ref,
                                 std::string strName,
                                 double      dDefaultValue,
                                 std::string strDescription) : CommandLineParameter(strName, strDescription),
                                                               m_dValue(ref)
      {
        m_dValue = dDefaultValue;
      }

      //  Identify the type of parameter.
      //===================================
      virtual ParameterType GetParameterType() const { return ParameterTypeDouble; }

      //  Number of values required.
      //==============================
      virtual int GetValueCount() const { return (1); }

      //  Print the usage information of this parameter.
      //=================================================
      virtual void PrintUsage() const 
      {
        printf("  %s <double> [%f] %s\n", m_strName.c_str(), m_dValue, 
                                          m_strDescription.c_str());
      }

      //  Set the value from a string.
      //===============================
      virtual void SetValue( int ix, std::string strValue) 
      {
        if (ix != 0) { _EXCEPTIONT("Invalid value index."); }
        m_dValue = atof(strValue.c_str());
      }

    public:
      double & m_dValue;   //  Parameter value.
};
//=====================================================================================


//=====================================================================================
//  Begin the definition of command line parameters.
//===================================================
#define BeginCommandLine() \
    { \
    bool _errorCommandLine = false; \
    std::vector<CommandLineParameter*> _vecParameters;

//  Define a new command line boolean parameter.
//===============================================
#define CommandLineBool(ref, name) \
    _vecParameters.push_back( \
      new CommandLineParameterBool(ref, name, ""));

//  Define a new command line boolean parameter with description.
//===============================================
#define CommandLineBoolD(ref, name, desc) \
    _vecParameters.push_back( \
      new CommandLineParameterBool(ref, name, desc));

//  Define a new command line string parameter.
//===============================================
#define CommandLineString(ref, name, value) \
    _vecParameters.push_back( \
      new CommandLineParameterString(ref, name, value, ""));

//  Define a new command line string parameter with description.
//===============================================
#define CommandLineStringD(ref, name, value, desc) \
    _vecParameters.push_back( \
      new CommandLineParameterString(ref, name, value, desc));

//  Define a new command line integer parameter.
//===============================================
#define CommandLineInt(ref, name, value) \
    _vecParameters.push_back( \
      new CommandLineParameterInt(ref, name, value, ""));

//  Define a new command line integer parameter with description.
//===============================================
#define CommandLineIntD(ref, name, value, desc) \
    _vecParameters.push_back( \
      new CommandLineParameterInt(ref, name, value, desc));

//  Define a new command line double parameter.
//===============================================
#define CommandLineDouble(ref, name, value) \
    _vecParameters.push_back( \
      new CommandLineParameterDouble(ref, name, value, ""));

//  Define a new command line double parameter with description.
//===============================================
#define CommandLineDoubleD(ref, name, value, desc) \
    _vecParameters.push_back( \
      new CommandLineParameterDouble(ref, name, value, desc));

//  Begin the loop for command line parameters.
//===============================================
#define ParseCommandLine(argc, argv) \
    for(int _command = 1; _command < argc; _command++) { \
      bool _found = false; \
      for(int _p = 0; _p < _vecParameters.size(); _p++) { \
        if (_vecParameters[_p]->m_strName == argv[_command]) { \
          _found = true; \
          _vecParameters[_p]->Activate(); \
          int _z; \
          if (_vecParameters[_p]->GetValueCount() >= argc - _command) { \
            printf("Error: Insufficient values for option %s\n", \
                    argv[_command]); \
            _errorCommandLine = true; \
            _command = argc; \
            break; \
          } \
          for (_z = 0; _z < _vecParameters[_p]->GetValueCount(); _z++) { \
            if ((_command + _z + 1 < argc) && \
                (strlen(argv[_command + _z + 1]) > 2) && \
                (argv[_command + _z + 1][0] == '-') && \
                (argv[_command + _z + 1][1] == '-') \
               ) break; \
            _command++; \
            _vecParameters[_p]->SetValue(_z, argv[_command]); \
          } \
          if ((_vecParameters[_p]->GetValueCount() >= 0) && \
              (_z != _vecParameters[_p]->GetValueCount()) \
             ) { \
            printf("Error: Insufficient values for option %s\n", \
                   argv[_command]); \
            _errorCommandLine = true; \
            _command = argc; \
            break; \
          } \
        } \
      } \
      if (!_found) { \
        printf("Error: Invalid parameter \"%s\"\n", argv[_command]); \
        _errorCommandLine = true; \
        break; \
      } \
    }

//  Print usage information.
//===============================================
#define PrintCommandLineUsage(argv) \
    if (_errorCommandLine) \
      printf("\nUsage: %s <Parameter List>\n", argv[0]); \
    printf("Input Parameters:\n"); \
    for (int _p = 0; _p < _vecParameters.size(); _p++) \
      _vecParameters[_p]->PrintUsage(); \
    if (_errorCommandLine) \
      exit(-1);

//  End the definition of command line parameters.
//===============================================
#define EndCommandLine(argv) \
    PrintCommandLineUsage(argv); \
    for (int _p = 0; _p < _vecParameters.size(); _p++) \
      delete _vecParameters[_p]; \
    }

#endif
