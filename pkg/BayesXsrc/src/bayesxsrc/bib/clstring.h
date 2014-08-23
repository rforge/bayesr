/* BayesX - Software for Bayesian Inference in
Structured Additive Regression Models.
Copyright (C) 2011  Christiane Belitz, Andreas Brezger,
Thomas Kneib, Stefan Lang, Nikolaus Umlauf

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */



#if !defined (CLSTRING_INCLUDED)
#define CLSTRING_INCLUDED

#include"../export_type.h"
#include<string.h>
#include<stdio.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<sys/stat.h>
#include<list>
#include<vector>


using std::vector;
using std::list;
using std::ostream;
using std::istream;

//------------------------------------------------------------------------------
//----------------------------- CLASS string -----------------------------------
//------------------------------------------------------------------------------

namespace ST
{

class __EXPORT_TYPE string
  {

  private:


  //------------------------- PRIVATE VARIABLES --------------------------------

  // stores the string

  char * str;

  // length of the string

  unsigned len;


  //------------------------- PRIVATE FUNCTIONS --------------------------------

  // FUNCTION: checkindex
  // TASK: checks, if access to charakter 'i' (of the calling string) is
  //       possible

  void checkindex(unsigned i) const
	 {
	 assert (len > 0);
	 assert(i < len);
	 }


  public:


  //------------------------- PUBLIC FUNCTIONS ---------------------------------

  // DEFAULT CONSTRUCTOR
  // ADDITIONAL INFORMATION:
  // - length = 0

  string();

  // CONSTRUCTOR
  // TASK: creates a string that consists of 'len' signs 'sign';

  string(const char & sign,const unsigned & len);

  // COPY CONSTRUCTORS

  string(const char * s);
  string(const string & s);
  string(const std::string & s);

  // DESTRUCTOR

  ~string()
	 {
	 delete [] str;
	 }

    // FUNCTION: strtochar
    // TASK: converts the calling string to char* and returns the result

    char * strtochar() const;

   // FUNCTION: to_bstr
   // TASK: converts the string to Std::string an returns the result

   std::string to_bstr(void) const;

  // OVERLOADED ASSIGNMENT OPERATORS

  const string & operator=(const string & s);
  const string & operator=(string & s);
  const string & operator=(const char * s);
  const string & operator=(const std::string & s);

  // OVERLOADED + OPERATORS

  friend string __EXPORT_TYPE operator+(const string & st1,const string & st2);
  friend string __EXPORT_TYPE operator+(const char * s,const string & st);

  // OVERLOADED OUTPUT OPERATOR

  friend ostream & __EXPORT_TYPE operator<<(ostream & c, const string & s)
	 {
	 c << s.str;
     return c;
	 }

  // OVERLOADED INPUT OPERATOR

  friend istream & __EXPORT_TYPE operator>>(istream & i, string & s);

  // FRIEND FUNCTION: getline
  // TASK: takes signs from the input stream 'i' until the delimeter
  //       is reached or number of signs is 'maxlen'.
  //       the resulting string is stored in 's'

  friend istream & __EXPORT_TYPE getline(istream & i,unsigned int maxlen,
									string & s, char delim);

  // FRIEND FUNCTION: getline
  // TASK: see function above
  // ADDITIONAL INFORMATION:
  // - maxlen = 256 Bytes

  friend istream & __EXPORT_TYPE getline(
  istream & i, string & s, char delim);

  // FRIEND FUNCTION: open
  // TASk:

//friend void open(std::ifstream & fin,string & s, int mode);

  friend void __EXPORT_TYPE open(
  std::ifstream & fin,string & s,int mode);

//MICRO  friend void open(ofstream & out,string & s,ios_base::openmode mode = ios_base::out);
  friend void __EXPORT_TYPE open(
  std::ofstream & out,string & s,int mode);


  // OVERLOADED [] OPERATOR
  // TASK: Zugriff auf ein einzelnes Zeichen
  // ADDITIONAL INFORMATION:  if i is invalid, the program terminates and an
  // errormessage occurs

  char & operator[] (int i);

  // OVERLOADED ASSIGNMENT OPERATORS

  friend int __EXPORT_TYPE operator==(const string & s1, const char * s2);

  friend int __EXPORT_TYPE operator==(const string & s1, const string & s2);

  friend int __EXPORT_TYPE operator!=(const string & s1, const char * s2);

  friend int __EXPORT_TYPE operator!=(string & s1, string & s2)
	 {
	 assert(s1.str != NULL);
	 assert(s2.str != NULL);
	 return strcmp(s1.str,s2.str) != 0;
	 }

  int __EXPORT_TYPE operator<(const string & s2) const
	  {
	  if (len < s2.len)
		 return 1;
	  else if (len > s2.len)
		 return 0;
	  else
		 {
       return strcmp(str,s2.str) < 0;
		 }
	  }

	friend int __EXPORT_TYPE operator>(string & s1,string & s2)
	  {
	  return s2 < s1;
     }


  // FUNCTION: length
  // TASK: returns the length of the string

  int length(void) const
	 {
	 return len;
	 }

  // Liefert die Länge des Teilsstrings der ausschließlich aus Zeichen besteht,
  // die in s enthalten sind

  int spn(const string & s) const
	 {
	 return strspn(str,s.str);
	 }

  int spn(char * s) const
	 {
	 return strspn(str,s);
	 }

	// FUNCTION: firstpos
	// TASK: returns the position in the calling string of the first appearence
	//       of 'sign', else -1

	int firstpos (char sign) const;

  // FUNCTION: helpfill
  // TASK: returns the string + the number of whitespaces that
  //       are necessary to complete the wanted width of the output column

  string helpfill(unsigned n);

  // FUNCTION: substr
  // TASK: returns a substring of the calling string
  //       copies 'nr' signs of the string starting at position 'pos'
  // POSSIBLE ERRORS:
  // - program terminates if  pos+nr > length of the string

  string substr(unsigned pos,unsigned nr) const;

  // FUNCTION: deletesign
  // TASK: deletetes the sign at position 'pos' of the calling string
  //       and returns the result
  // ADDITIONAL INFORMATION:
  // asserts if pos is a valid index (leads to an assertion failure,
  // abnormal program termination)

  string deletesign(unsigned pos) const;

  // FUNCTION: deleteallsigns
  // TASK: deletes all 'sign' signs of the calling string and returns the
  //       result
  // EXAMPLE:
  // string test = "today"
  // test.deleteallsigns('d') returns the string "toay"

  string deleteallsigns(char sign) const;

  // FUNCTION: replaceallsigns
  // TASK: replaces all 'oldsign' signs of the calling string with 'newsign' signs
  //       and returns the result
  // EXAMPLE:
  // string test = "today"
  // test.replaceallsigns('d','k') returns the string "tokay"

  string replaceallsigns(char oldsign, char newsign) const;

  //FUNCTION: insert_string_num
  //TASK: inserts the string 'str' at the position 'pos' into a string
  //EXAMPLE:
  // string test = "today"
  // test.insert_string(2, "hi") returns the string "tohiday"

  string insert_string_num(unsigned pos, string & str) const;

  //FUNCTION: insert_string_char
  //TASK: inserts the string 'str' instead of each character 'p' into a string
  //EXAMPLE:
  // string test = "today"
  // test.insert_string('o', "hi") returns the string "thiday"

  string insert_string_char(char p, string & str) const;

  //FUNCTION: insert_after_string
  //TASK: inserts the string 's1' after the first appearance of the string 's2'
  //NOTE: returns the calling string if there is no appearance of 's2'
  //EXAMPLE:
  // string test = "today"
  // test.insert_after_string("hi", "od") returns the string "todhiay"

  string insert_after_string(string s1, string s2) const;

  //FUNCTION: insert_after_all_string
  //TASK: inserts the string 's1' after every appearance of the string 's2'
  //NOTE: returns the calling string if there is no appearance of 's2'
  //EXAMPLE:
  // string test = "today"
  // test.insert_after_string("hi", "od") returns the string "todhiay"

  string insert_after_all_string(string s1, string s2) const;

  // FUNCTION: eatwhitespace
  // TASK: returns a string, where (possible) leading and ... whitespace signs
  //       of the calling string are deleted
  // EXAMPLE:
  // string test = "  today ";
  // string result = test.eatwhitespace();
  // result evaluates to "today"

  string eatwhitespace(void) const;

  // FUNCTION: eatallwhitespace
  // TASK: deletes all whitespacesigns of the calling string and returns the
  //       result
  // EXAMPLE:
  // string test = " hello world"
  // test.eatallwhitespace() returns the string "helloworld"

  string eatallwhitespace(void) const;

  string eatallcarriagereturns(void) const;

  // FUNCTION: closingbracketpos
  // TASK: returns the position of the closing bracket ) in the calling string,
  //       if 'bracketpos' is the position of an opening bracket
  //       returns -1 if the closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - asserts 0 <= bracketpos < length  of the string
  // - asserts if position 'bracketpos' contains an opening bracket '('
  // EXAMPLE:
  // string test = "sin(x+y)+4"
  // test.closingbracketpos(3) returns 7
  // string test = "sin(4+x"
  // test.closingbracketpos(3) returns -1 (closing bracket is missing)

  int closingbracketpos(const unsigned bracketpos) const;

  // FUNCTION: closingbracketpos
  // TASK: returns the position of the closing bracket ] in the calling string,
  //       if 'bracketpos' is the position of an opening bracket
  //       returns -1 if the closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - asserts 0 <= bracketpos < length  of the string
  // - asserts if position 'bracketpos' contains an opening bracket '('
  // EXAMPLE:
  // string test = "sin[x+y]+4"
  // test.closingbracketpos(3) returns 7
  // string test = "sin[4+x"
  // test.closingbracketpos(3) returns -1 (closing bracket is missing)

  int closingbracketpos2(const unsigned bracketpos) const;

  // FUNCTION: lowestprecedencepos
  // TASK:

  int lowestprecedencepos(string & sign) const;


  // FUNCTION: isfunction
  // TASK: returns 1, if the calling string is a function, i.e. of type
  //       functionname(argument)
  //       in 'functionname' the name of the function will be stored,
  //       in 'argument' the argument of the function will be stored
  //       returns 0, if calling string is not a function
  //        returns -1 if closing bracket is missing
  // ADDITIONAL INFORMATION:
  // - whitespace signs should be deleted, before calling the function
  // EXAMPLE:
  // sin(x+z-3).isfunction(functionname,argument) will return 1 and
  // functionname = "sin" and argument = "x+z-3"

  int isfunction(string & functionname,string & argument) const;

  int issubscribing(string & varname, string & argument) const;

  // FUNCTION: isexistingfile
  // TASK: checks, if calling string contains a path to an existing file,
  //       that can be opened for reading
  //       returns 0 if path is valid
  //               1 if path is invalid (error)

  int isexistingfile(void) const;

  // FUNCTION: isvalidfile
  // TASK: checks, if calling string contains a path to a file, that can be
  //       opened for writing
  //       returns: 1 if file can not be opened for writing
  //                0 if file can be opened for writing and is not existing
  //                -1 if file can be opened for writing but is already existing

  int isvalidfile(void) const;


  // FUNCTION: isvarname
  // TASK: checks, if calling string is a  valid variable name
  //       returns: 1 = error, invalid variablename
  //                0 = no error
  // VALID VARIABLE NAME:
  // - first sign is a literal
  // (abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_)
  // - following signs are either literals or numbers (0123456789)

  int isvarname() const;



  int removefile(void) const
	 {
	 return remove(str);
	 }


  int isint(void) const;

  // FUNCTIONS: strtolong,strtochar,strtouchar,strtodouble
  // TASK: functions, that convert strings into numbers
  //       Converts the calling string into the long value
  //       (char value, unsigned char value, double value) "value"
  //       returns:  0 if no error occured
  //                 1 = if an error occured (value will not be changed)

  int strtolong(long & value) const;

  int strtochar(char & value) const;

  int strtouchar(unsigned char & value) const;

  int strtodouble(double & value) const;

  // FUNCTION: checksign
  // TASK: checks, if the sign 'signs' is a member of the calling string
  //       returns the position of the first sign in the calling string,
  //       that is equal to 'signs'
  //       returns -1 if 'sign' is not a member of the calling string

  int checksign(const char sign) const;

  // FUNCTION: isinlist
  // TASK: returns -1, if calling string is not a meber of 'stringlist'
  //                position of the string in stringlist otherwise

  int isinlist(const vector<string> & stringlist) const;

  // FUNCTION: getFirstToken
  // TASK: returns the first token of the calling string
  //       if length = 0 an empty string will be returned
  //       delimeters are stored in 'parsingsigns'

  string getFirstToken(const string & parsingsigns) const;

  // FUNCTION: strtoken
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingsigns'

  vector<string> strtoken(const string & parsingsigns,bool inclsigns = true) const;

  int strtoken_quot(vector<string> & hilfe,const string & parsingsigns,
  bool inclsigns = true) const;

  // FUNCTION: strtoken2
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingsigns'
  // ADDITIONAL INFORMATION:
  // Difference to strtoken: signs in brackets are ignored

  vector<string> strtoken2(const string & parsingsigns,bool & bracketmiss) const;

  vector<string> strtoken2_quot(const string & parsingsigns,
                                bool & bracketmiss, bool & quotmiss) const;

  // FUNCTION: strtoken
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingtokens'

  vector<string> strtoken(const vector<string> & parsingtokens) const;

  // FUNCTION: strtoken2
  // TASK: returns the tokens of the string as a vector of strings
  //       delimeters are stored in 'parsingtokens'
  // ADDITIONAL INFORMATION:
  // Difference to strtoken: signs in brackets are ignored
  // If closing brackets are missing, bracketmiss = true, else false

  vector<string> strtoken2(const vector<string> & parsingtokens,
                           bool & bracketmiss) const;


  vector<string> strtoken2_quot(const vector<string> & parsingtokens,
                                      bool & bracketmiss,bool & quotmiss) const;


  list<string> strtokenlist(const string & parsingsigns,bool  inclsigns = true) const;

  // FUNCTION: endswith
  // TASK: checks if the calling string ends with 'c'

  bool endswith(const char * c) const;


  };  // end: class string



string __EXPORT_TYPE outresults(const unsigned & l,const string & name,
                  const double & mean,
                  const double & std, const double & qu10, const double & qu50,
                  const double & qu90);

string __EXPORT_TYPE make_latextable(vector<string> & v);

//------------ forward declarations of friends ---------------------------------
istream & getline(istream & i,string & s, char delim = '\n');
istream & getline(istream & i,unsigned int maxlen,string & s,char delim = '\n');

void __EXPORT_TYPE open (std::ifstream & fin, string & s, int mode = std::ios::in);
void __EXPORT_TYPE open(std::ofstream & out,string & s,int mode = std::ios::out);

//------------ functions, that convert a number into a string ------------------

// FUNCTION: inttostring
// TASK: converts the integer number 'value' into a string

string __EXPORT_TYPE inttostring(int value);

// FUNCTION: doubletostring
// TASK: converts the double number 'value' into a string

string __EXPORT_TYPE doubletostring(double value,int dec=15);


//------------------------ type definitons -------------------------------------

typedef vector<string> __EXPORT_TYPE stringvec;
typedef list<string> __EXPORT_TYPE stringlist;

/*
template<class T>
void pr(const vector<T> & v)
  {
  if (! v.empty())
	 for (int i = 0;i < v.size();i++)
		cout << v[i] << "\n";
  cout << "\n";
  }

template<class T>
void pr(list<T> & v)
  {
  if (!v.empty())
	 {
	 list<T>::iterator i;
	 for(i=v.begin();i != v.end();i++)
		cout << (*i) << endl;
	 }
  }
*/


}  // end: namespace ST

#endif


