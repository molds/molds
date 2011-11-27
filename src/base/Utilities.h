#ifndef INCLUDED_UTILITIES
#define INCLUDED_UTILITIES
namespace MolDS_base{

/* 日時を "hoge/hoge/hoge 00:00:00" にフォーマット 
 *      引数 : buf...文字列の格納場所、st...参照するtm構造体へのポインタ
 *           戻り値 : 格納した文字列へのポインタ        
 *
 *  http://www.sasaraan.net/program/cpp/cpp_time.html
 *           */
char *fmttm(char *buf, struct tm *st);

// trim the string
std::string TrimString(const std::string str);

// number to string
// ex. Num2String(23,5) = "00023";
std::string Num2String(int number, int digit);

}
#endif

