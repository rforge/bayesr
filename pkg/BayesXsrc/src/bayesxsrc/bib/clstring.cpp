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
#include "clstring.h"
#endif

#include <stdlib.h>
#include <climits>
#include <sstream>
#include <string>

#include <iomanip>
#include <iostream>

using std::ifstream;
using std::ofstream;
using std::ios;


//------------------------------------------------------------------------------
//-------------- CLASS STRING: Implementation of member functions --------------
//------------------------------------------------------------------------------

namespace ST
{

string::string()
  {
  len = 0;
  str = new char[1];
  strcpy(str,"");
  }

string::string(const char & sign,const unsigned & l)
  {
  len = l;
  str = new char[len+1];
  for (unsigned i=0;i<len;i++)
	 str[i] = sign;
  str[len] = '\0';
  }


string::string(const string & s)
  {
  len = s.len;
  str = new char[len+1];
  strcpy(str,s.str);
  }


string::string(const char * s)
  {
  len = strlen(s);
  str = new char[len+1];
  strcpy(str,s);
  }

string::string(const std::string & s)
  {
  len = s.length();
  str = new char[len+1];
  strcpy(str,s.c_str());
  }


std::string string::to_bstr(void) const
  {
  return std::string(str);
  }

void open(ifstream & fin,string & s,int mode)
  {
// GNU:
  fin.open(s.strtochar(),ios::in);
//  fin.open(s.strtochar(),mode);
  }


void open(ofstream & out, string & s,int mode)
  {
// GNU:
  out.open(s.strtochar(),ios::out);
//  out.open(s.strtochar(),mode);
  }

const string & string::operator=(const string & s)
  {
  if (this == &s)
	 return *this;
  delete [] str;
  len = s.len;
  str = new char[len+1];
  strcpy(str,s.str);
  return *this;
  }


const string & string::operator=(string & s)
  {
  if (this == &s)
	 return *this;
  delete [] str;
  len = s.len;
  str = new char[len+1];
  strcpy(str,s.str);
  return *this;
  }


const string & string::operator=(const char * s)
  {
  delete [] str;
  len = strlen(s);
  str = new char[len+1];
  strcpy(str,s);
  return *this;
  }


const string & string::operator=(const std::string & s)
  {
  delete[] str;
  len = s.length();
  str = new char[len+1];
  strcpy(str,s.c_str());
  return *this;
  }


int operator==(const string & s1, const char * s2)
  {
  assert(s1.str != NULL);
  assert(s2 != NULL);
  return strcmp(s1.str,s2) == 0;
  }


int operator==(const string & s1, const string & s2)
  {
  assert(s1.str != NULL);
  assert(s2.str != NULL);
  return strcmp(s1.str,s2.str) == 0;
  }


int operator!=(const string & s1, const char * s2)
  {
  assert(s1.str != NULL);
  assert(s2 != NULL);
  return strcmp(s1.str,s2) != 0;
  }


string operator+(const string & st1,const string & st2)
  {
  char * h = new char[st1.len+st2.len+1];
  strcpy(h,st1.str);
  strcpy(h+st1.len,st2.str);
  string sum(h);
  delete [] h;
  return sum;
  }


string operator+(const char * s,const string & st)
  {
  string h(s);
  return h+st;
  }


istream & operator>>(istream & i, string & s)
  {
  char buffer[256];
  i >> buffer;
  s = buffer;
  return i;
  }


istream & getline(istream & i,unsigned int maxlen,string & s,char delim)
  {
  char * buffer = new char[maxlen];
  i.getline(buffer,maxlen,delim);
  s = buffer;
  delete [] buffer;
  return i;
  }


istream & getline(istream & i,string & s,char delim)
  {
  return getline(i,4096,s,delim);
  }


char & string::operator[](int i)
  {
  checkindex(i);
  return str[i];
  }


string string::helpfill(unsigned n)
  {
  unsigned d;
  ST::string wert;
  if (n>=len)
    {
    d = n - len;
    wert = substr(0, len);
    }
  else
    {
    wert = substr(0, n-2) + "~";
    d = 1;
    }
  ST::string platz = string(' ',d);
  ST::string out = platz + wert;
  return out;
  }


string string::substr(unsigned pos, unsigned nr) const
  {
  assert(pos+nr <= len);
  assert(nr > 0);

  char * help = new char[nr+1];
  strncpy(help,str+pos,nr);
  help[nr] = '\0';
  string ret(help);
  delete [] help;
  return ret;
  }


int string::firstpos (char sign) const
  {
  int pos = -1;
  unsigned i = 0;
  while ( (i < len) && (pos == -1) )
	 {
	 if (str[i] == sign)
		pos = i;
	 i++;
	 }
  return pos;
  }


string string::deletesign(unsigned pos) const
  {
  checkindex(pos);
  if (pos == 0)
    {
    if (len > 1)
      return substr(1,len-1);
    else
      return "";
    }
  else if (pos == len-1)
	 return substr(0,len-1);
  else
    return substr(0,pos) + substr(pos+1,len-pos-1);
  }


string string::deleteallsigns(char sign) const
  {
  string result = *this;
  int i=0;
  while ( i < result.length() )
	 {
	 if (result[i] == sign)
		{
		result = result.deletesign(i);
		}
	 else
		i++;
	 }
  return result;
  }



string string::replaceallsigns(char oldsign, char newsign) const
  {
  int i=0;
  string result = *this;
  while (i< result.length())
	 {
	 if (result[i] == oldsign)
		result[i] = newsign;
	 i++;
	 }
  return result;
  }



string string::insert_string_num(unsigned pos, string & str) const
  {
  string s = *this;
  assert(pos<s.length());
  string s1 = s.substr(0, pos);
  string s2 = s.substr(pos, s.length()-pos);
  string result = s1 + str + s2;
  return result;
  }



string string::insert_string_char(char p, string & str) const
  {
  string s = *this;
  string result = s;
  unsigned l = str.length();
  unsigned k = 0;
  for (unsigned i=0;i<s.length()-1;i++)
      {
      char z = s[i];
      if (z == p)
         {
         string s1 = result.substr(0, i+k*l-k);
         string s2 = result.substr(i+1+k*l-k, result.length()-(i+1+k*l-k));
         result = s1 + str + s2;
         k = k + 1;
         }
      }
  return result;
  }

string string::insert_after_string(string s1, string s2) const
  {
  string s = *this;
  string shelp;
  unsigned n2 = s2.length();
  unsigned n = s.length()-n2;
  unsigned i;
  for(i=0; i<n+1; i++)
    {
    shelp=s.substr(i,n2);
    if(shelp==s2)
      {
      shelp = s.substr(0,i+n2)+s1;
      if(i<n-1)
        {
        shelp = shelp+s.substr(i+n2,n-i);
        }
      return shelp;
      }
    }
  return s;
  }

string string::insert_after_all_string(string s1, string s2) const
  {
  string s = *this;
  string shelp;
  string result = " ";

  unsigned n2 = s2.length();
  unsigned n = s.length()-n2;
  unsigned i=0,k=0;
  unsigned j=0;

  for(i=0; i<n+1; i++)
    {
    shelp=s.substr(i,n2);
    if(shelp==s2)
      {
      if(j==0)
        {
        result=s.substr(0,i+n2)+s1;
        j++;
        }
      else
        {
        result = result+s.substr(k+n2,i-k)+s1;
        }
      k=i;
      }
    }
  if(j==1)
    {
    if(k<n-1)
      {
      result = result + s.substr(k+n2,n-k);
      }
    return result;
    }
  return s;
  }

string string::eatwhitespace(void) const
  {

  int beginpos;
  int endpos;

  int i=0;

  while ( (i<len) && (str[i] == ' ') )
	 i++;
  if (i == len)
	 return "";
  beginpos = i;

  i=len-1;
  while ( (i >= 0) && (str[i] == ' ') )
	 i--;
  endpos=i;

  return substr(beginpos,endpos-beginpos+1);

  }

string string::eatallwhitespace(void) const
  {
/*  string res = *this;
  res.deleteallsigns(' ');
  res.deleteallsigns('\n');
  res.deleteallsigns('\r');
  return res;*/
  return deleteallsigns(' ');
  }

string string::eatallcarriagereturns(void) const
  {
  return deleteallsigns('\r');
  }


int string::closingbracketpos(const unsigned bracketpos) const
  {
  assert (bracketpos < len);
  assert (str[bracketpos] == '(');

  unsigned i=bracketpos+1;
  unsigned nropen = 1;
  while ( (i < len) && (nropen > 0) )
	 {
	 if (str[i] == '(')
		nropen++;
	 else if (str[i] == ')')
		nropen--;
	 i++;
	 }
  if (nropen == 0)
	 return i-1;
  else
	 return -1;
  }

int string::closingbracketpos2(const unsigned bracketpos) const
  {
  assert (bracketpos < len);
  assert (str[bracketpos] == '[');

  unsigned i=bracketpos+1;
  unsigned nropen = 1;
  while ( (i < len) && (nropen > 0) )
	 {
	 if (str[i] == '[')
		nropen++;
	 else if (str[i] == ']')
		nropen--;
	 i++;
	 }
  if (nropen == 0)
	 return i-1;
  else
	 return -1;
  }



int string::lowestprecedencepos(string & sign) const
  {
  int i = 0;
  int pos = -1;
  int minp = 6;
  while (i < len)
	 {
	 if (str[i] == '(')
		{
		i = closingbracketpos(i);
		if (i == -1)
		  return -2;
		i++;
		}
     else if (str[i] == '[')
       {

		i = closingbracketpos2(i);
		if (i == -1)
		  return -2;
		i++;

       }
	 else if ( (str[i] == '+') || (str[i] == '-') )
		{
		if (minp >= 3)
		  {
		  minp = 3;
		  pos = i;
		  if (str[i] == '+')
			 sign = "+";
		  else
			 sign = "-";
		  }
		i++;
		}
	 else if ( (str[i] == '*') || (str[i] == '/') )
		{
		if (minp >= 4)
		  {
		  minp = 4;
		  pos = i;
		  if (str[i] == '*')
			 sign = "*";
		  else
			 sign = "/";
		  }
		i++;
		}
	 else if ( (str[i] == '^') )
		{
		if (minp >= 5)
		  {
		  minp = 5;
		  pos = i;
		  sign = "^";
		  }
		i++;
		}
	 else if (str[i] == '=')
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = "=";
		  }
		i++;
		}
	 else if ( (str[i] == '>') && (i+1 < len) && (str[i+1] == '=') )
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = ">=";
		  }
		i=i+2;
		}
	 else if ( (str[i] == '<') && (i+1 < len) && (str[i+1] == '=') )
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = "<=";
		  }
		i=i+2;
		}
	 else if ( (str[i] == '!') && (i+1 < len) && (str[i+1] == '=') )
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = "!=";
		  }
		i=i+2;
		}
	 else if (str[i] == '>')
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = ">";
		  }
		i++;
		}
	 else if (str[i] == '<')
		{
		if (minp >= 2)
		  {
		  minp = 2;
		  pos = i;
		  sign = "<";
		  }
		i++;
		}
	 else if (str[i] == '&')
		{
		if (minp >= 1)
		  {
		  minp = 1;
		  pos = i;
		  sign = "&";
		  }
		i++;
		}
	 else if (str[i] == '|')
		{
		if (minp >= 1)
		  {
		  minp = 1;
		  pos = i;
		  sign = "|";
		  }
		i++;
		}
	 else
		i++;
	 }
  return pos;
  }



int string::isfunction(string & functionname,string & argument) const
  {
  int startbr = firstpos('(');
  if (startbr > 0)
    {
    int endbr = closingbracketpos(startbr);
    if (endbr == -1)
      return -1;
    if (endbr == len-1)
      {
      functionname = substr(0,startbr);
      if (len -startbr-2 > 0)
        argument = substr(startbr+1,len-startbr-2);
      else
        argument = "";
      return 1;
      }
    else
      return 0;
    }
  else
	 return 0;
  }




int string::issubscribing(string & varname,string & argument) const
  {
  int startbr = firstpos('[');
  if ( (startbr > 0) && (startbr != -1) )
	 {
	 int endbr = closingbracketpos2(startbr);
	 if (endbr == len-1)
		{
		varname = substr(0,startbr);
		argument = substr(startbr+1,len-startbr-2);
		return 1;
		}
	 else
		return -1;
	 }
  else
	 return -1;
  }




int string::isexistingfile(void) const
  {
  ifstream fin(str,ios::in);
  if (fin.fail() != 0)
	 return 1;
  else
	 return 0;
  }


int string::isvalidfile(void) const
  {
  struct stat statbuf;
  int existing = stat(str,&statbuf);
  if (existing == 0)
	 {
	 ofstream fout(str,ios::app);
	 if (fout.fail())
		return 1;
	 else
		return -1;
	 }
  else
	 {
	 ofstream fout(str);
	 if (fout.fail())
		{
		fout.close();
		remove(str);
		return 1;
		}
	 else
		{
		fout.close();
		remove(str);
		return 0;
		}
	 }
  }


int string::isvarname() const
  {
  if (len > 0)
	 {
	 string valid =
	 "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_0123456789";
	 unsigned anz = spn(valid);                       // number of valid signs
	 if (anz < len)
		return 1;
	 else
		{
		string numbers = "0123456789";
		string firstsign = substr(0,1);
		if (firstsign.spn(numbers) != 0)
		  return 1;
		else
		  return 0;
		}
	 }
  else
	 return 1;
  }


int string::isint(void) const
  {
  if ((len == 0))
	 return 0;
  else
	 {
	 int h = 1;
	 int i;
	 if ((str[0] == '-')  || (str[0] == '+'))
		i = 1;
	 else
		i = 0;
	 string firstsign = "123456789";
	 string nextsigns = "0123456789";
	 if (firstsign.checksign(str[i]) == -1)
		h = 0;
	 i++;
	 while ((i < len) && (h == 1))
		{
		if (nextsigns.checksign(str[i]) == -1)
		  h = 0;
		i++;
		}
	 return h;

	 }
  }


char * string::strtochar() const
  {
  char * h = new char[len+1];
  strcpy(h,str);
  return h;
  }


int string::strtolong(long & value) const
  {
  if (len > 0)
	 {
	 char * sentinel;
	 long h = strtol(str, & sentinel,10);
	 if (sentinel != str+len)
		return 1;
	 else
		{
		value = h;
		return 0;
		}
	 }
  else
	 return 1;
  }


int string::strtochar(char & value) const
  {
  long h;
  if (strtolong(h) == 1)
	 return 1;
  else if ((h < CHAR_MIN) || (h > CHAR_MAX))
	 return 1;
  else
	 {
	 value = h;
	 return 0;
	 }
  }


int string::strtouchar(unsigned char & value) const
  {
  long h;
  if (strtolong(h) == 1)
	 return 1;
  else if ((h < 0) || (h > UCHAR_MAX))
	 return 1;
  else
	 {
	 value = h;
	 return 0;
	 }
  }


int string::strtodouble(double & value) const
  {
  if (len > 0)
	 {
	 char * sentinel;
	 double h = strtod(str, & sentinel);
	 if (sentinel != str+len)
		return 1;
	 else
		{
		value = h;
		return 0;
		}
	 }
  else
	 return 1;
  }


string inttostring(int value)
  {
  char h[20];
// GNU:
  sprintf(h,"%d",value);
//  itoa(value,h,10);
  return string(h);
  }

string doubletostring(double value,int dec)
  {
  if ((dec > 19) || (dec < 1))
     {
	 dec = 15;
     }
/*  char h[20];
  gcvt(value,dec,h);
  string s = h;
  return s;*/
  std::ostringstream oss;
  oss << std::setprecision(dec) << value;
  return oss.str();
//  return string(h);
  }

int string::checksign(const char sign) const
  {
  unsigned i=0;
  int isfrom = -1;
  while ((i < len) && (isfrom == -1))
	 {
	 if (str[i] == sign)
		isfrom = i;
	 i++;
	 }
  return isfrom;
  }


int string::isinlist(const vector<string> & stringlist) const
  {
  unsigned i = 0;
  int isinl = -1;
  while ((i < stringlist.size()) && (isinl == -1))
	 {
	 if (*this == stringlist[i])
		isinl = i;
	 i++;
	 }
  return isinl;
  }


string string::getFirstToken(const string & parsingsigns) const
  {
  if (len > 0)
	 {
	 unsigned i;
	 unsigned j=0;
	 while (str[j] == ' ')
		j++;
	 i = j;
	 while ( (i < len) && (parsingsigns.checksign(str[i]) == -1))
		i++;
     if (i > j)
	   return substr(j,i-j);
     else
       return string();
	 }
  else
	 return string();
  }


vector<string> string::strtoken(const string & parsingsigns,bool inclsigns) const
  {
  vector<string>  hilfe;


  if (len > 0)
	 {
	 int i=0;                             // looping variable
	 while (i < len)
		{
		if (parsingsigns.checksign(str[i]) != -1)   // str[i] is a parsingsign
		  if (str[i] == ' ')
			 while ( (i < len) && (str[i] == ' ') )
				i++;
		  else
			 {
             if (inclsigns)
			   hilfe.push_back(substr(i,1));
			 i++;
			 }
		else                                       // str[i] is not a parsingsign
		  {
		  int anf = i;
		  while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
			 i++;
		  hilfe.push_back(substr(anf,i-anf));
		  }
		}
	 }
  return hilfe;
  }


int string::strtoken_quot(vector<string> & hilfe, const string & parsingsigns,
                                     bool inclsigns) const
  {

  bool ok = true;

  if (len > 0)
	 {
	 int i=0;                             // looping variable
	 while (i < len)
		{
		if (parsingsigns.checksign(str[i]) != -1)   // str[i] is a parsingsign
		  if (str[i] == ' ')
			 while ( (i < len) && (str[i] == ' ') )
				i++;
		  else
			 {
             if (inclsigns)
			   hilfe.push_back(substr(i,1));
			 i++;
			 }
		else                                       // str[i] is not a parsingsign
		  {
		  int anf = i;

          if (str[i] != '"')
            {
		    while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
              i++;
		    hilfe.push_back(substr(anf,i-anf));
            }
          else
            {
            i++;
            ok = false;
            while ( ( i< len) && (str[i] != '"') )
              i++;

            if ( ( i < len) && (str[i] == '"') )
              ok = true;

            i++;
            if (i-anf-2 > 0)
              hilfe.push_back(substr(anf+1,i-anf-2));
            else
              hilfe.push_back("");
            }

		  }
		}
	 }

  return ok;

  }


vector<string> string::strtoken2(const string & parsingsigns,
                                      bool & bracketmiss) const
  {
  vector<string>  hilfe;

  bracketmiss = false;

  if (len > 0)
	 {
	 int i=0;                             // looping variable
	 int anf = 0;
	 while (i < len)
		{
		if (str[i] == '(')
		  {
		  i = closingbracketpos(i);
		  if (i == -1)
            {
            bracketmiss = true;
            return vector<string>();
            }
		  else
			 i++;
		  }
		else if (parsingsigns.checksign(str[i]) != -1)
		  {
          if (i-anf > 0)
		    hilfe.push_back(substr(anf,i-anf));
		  i++;
		  anf = i;
		  }
		else
		  i++;
		}  // end: while (i <len)
	 if (anf < len)
		hilfe.push_back(substr(anf,len-anf));
	 }

  return hilfe;

  }



vector<string> string::strtoken2_quot(const string & parsingsigns,
                                      bool & bracketmiss,bool & quotmiss) const
  {
  vector<string>  hilfe;

  bracketmiss = false;
  quotmiss = false;

  if (len > 0)
	 {
	 int i=0;                             // looping variable
	 int anf = 0;
	 while (i < len)
		{
		if (str[i] == '(')
		  {
		  i = closingbracketpos(i);
		  if (i == -1)
            {
            bracketmiss = true;
            return vector<string>();
            }
		  else
			 i++;
		  }
        else if (str[i] == '"')
          {
          i++;
          while ( ( i< len) && (str[i] != '"') )
            i++;

          if ( ( i < len) && (str[i] == '"') )
            {
            }
          else
            quotmiss = true;

          i++;
          }
		else if (parsingsigns.checksign(str[i]) != -1)
		  {
          if (i-anf > 0)
		    hilfe.push_back(substr(anf,i-anf));
		  i++;
		  anf = i;
		  }
		else
		  i++;
		}  // end: while (i <len)
	 if (anf < len)
		hilfe.push_back(substr(anf,len-anf));
	 }

  return hilfe;

  }


vector<string> string::strtoken(const vector<string> & parsingtoken) const
  {

  vector<string> hilfe;

  vector<string> token = strtoken(" ");

  unsigned i = 0;
  while (i < token.size())
	 {
	 string h = token[i];
	 i++;
	 if (h.isinlist(parsingtoken) == -1)
		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
		  {
		  h = h + " " + token[i];
		  i++;
		  }
	 hilfe.push_back(h);
	 }
  return hilfe;
  }


vector<string> string::strtoken2_quot(const vector<string> & parsingtoken,
                                      bool & bracketmiss, bool & quotmiss) const
  {
  vector<string> hilfe;

  vector<string> token = strtoken2_quot(" ",bracketmiss,quotmiss);

  if ((!bracketmiss) && (!quotmiss))
   {
   unsigned i = 0;
   while (i < token.size())
     {
	 string h = token[i];
	 i++;
	 if (h.isinlist(parsingtoken) == -1)
		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
		  {
		  h = h + " " + token[i];
		  i++;
		  }
	 hilfe.push_back(h);
	 }
   }

  return hilfe;

  }


vector<string> string::strtoken2(const vector<string> & parsingtoken,
                                 bool & bracketmiss) const
  {
  vector<string> hilfe;

  vector<string> token = strtoken2(" ",bracketmiss);

  if (!bracketmiss)
  {
  unsigned i = 0;
  while (i < token.size())
	 {
	 string h = token[i];
	 i++;
	 if (h.isinlist(parsingtoken) == -1)
		while ((i < token.size()) && (token[i].isinlist(parsingtoken) == -1))
		  {
		  h = h + " " + token[i];
		  i++;
		  }
	 hilfe.push_back(h);
	 }
   }

  return hilfe;

  }


list<string> string::strtokenlist(const string & parsingsigns,bool inclsigns) const
  {
  list<string> help;                     // will be returned at the
													  // of the function
  if (len > 0)
	 {
	 int i=0;                           // looping variable
     while (i < len)
		{
		if (parsingsigns.checksign(str[i]) != -1)
		  if (str[i] == ' ')
			 while ( (i < len) && (str[i] == ' ') )
				i++;
		  else
			 {
             if (inclsigns)
			   help.push_back(substr(i,1));
			 i++;
			 }
		else
		  {
		  int anf = i;
		  while ( (i < len) && (parsingsigns.checksign(str[i]) == -1) )
            i++;

		  help.push_back(substr(anf,i-anf));
		  }
		}  // end: while (i < text.length())
	 }  // end: if (text.length() > 0)
  return help;
  }


bool string::endswith(const char * c) const
    {
    int lenc = strlen(c);
//    bool endwith = true;
    for(int i=0;i<lenc;i++)
      if(c[lenc-1-i] != str[len-1-i])
        return false;
    return true;
    }


string outresults(const unsigned & l,const string & name, const double & mean,
                  const double & std, const double & qu10, const double & qu50,
                  const double & qu90)
  {


  string means = ST::doubletostring(mean,6);
  means = means + ST::string(' ',15-means.length());

  string stds = ST::doubletostring(std,6);
  stds = stds + ST::string(' ',15-stds.length());

  string qu10s = ST::doubletostring(qu10,6);
  qu10s = qu10s + ST::string(' ',15-qu10s.length());

  string qu50s = ST::doubletostring(qu50,6);
  qu50s = qu50s + ST::string(' ',15-qu50s.length());

  string qu90s = ST::doubletostring(qu90,6);
  qu90s = qu90s + ST::string(' ',15-qu90s.length());

  string ls(' ',l);

  string result = "    " + name + ls + means + stds + qu10s + qu50s + qu90s;

  return result;

  }

string make_latextable(vector<string> & v)
  {
  string result="";
  unsigned i;

  for (i=0;i<v.size();i++)
    {
    if (i==0)
      result = v[i];
    else
      result = result + " & " + v[i];
    }

  result = result + "\\\\";

  return result;
  }


} // end: namespace ST
