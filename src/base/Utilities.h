#ifndef INCLUDED_UTILITIES
#define INCLUDED_UTILITIES
namespace MolDS_base{
// string for today.
std::string GetDateString();
// trim the string
std::string TrimString(const std::string str);
// number to string
// ex. Num2String(23,5) = "00023";
std::string Num2String(int number, int digit);

}
#endif

