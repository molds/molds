#include<stdio.h>
#include<stdlib.h>
#include<sstream>
#include<string>
#include<math.h>
#include"Utilities.h"
using namespace std;

namespace MolDS_base{
// string of today
string GetDateString(){
   time_t current;
   struct tm *local;
   char  wday_name[][10] = {"Sun.", "Mon.", "Thu.", "Wed.", "Thu.", "Fri.", "Sat."};
   time(&current);
   local = localtime(&current);
   stringstream ss;
   ss << local->tm_year + 1900 
      << "/" 
      << local->tm_mon + 1 
      << "/" 
      << local->tm_mday 
      << "(" 
      << wday_name[local->tm_wday] 
      << ") ";
   ss << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec;
   return ss.str();
}

// trim the string
string TrimString(const string str){
   int nStart = 0;
   int nEnd = str.length() - 1;
   // left trim 
   for(int n = 0; n < str.length(); n++ ){
      if( str.data()[n] != ' ' ){
         nStart = n;
         break;
      }
   }
   // right trim 
   for(int n = str.length() - 1; n >= 0; n-- ){
      if( str.data()[n] != ' ' ){
         nEnd = n;
         break;
      }
   }
   return(str.substr( nStart, nEnd - nStart + 1 ));
}

string Num2String(int number, int digit){
   stringstream ss;
   int numberDigit = (int)(log10((double)number)) + 1;
   for(int i=0; i<digit-numberDigit; i++){
      ss << "0";
   }
   ss << number;
   return ss.str();
}
}

