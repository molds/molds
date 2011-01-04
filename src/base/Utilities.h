#ifndef INCLUDED_UTILITIES
#define INCLUDED_UTILITIES

using namespace std;

namespace MolDS_base{

/* 日時を "hoge/hoge/hoge 00:00:00" にフォーマット 
 *      引数 : buf...文字列の格納場所、st...参照するtm構造体へのポインタ
 *           戻り値 : 格納した文字列へのポインタ        
 *
 *  http://www.sasaraan.net/program/cpp/cpp_time.html
 *           */
char *fmttm(char *buf, struct tm *st) 
{

    sprintf(buf, "%d/%d/%d %02d:%02d:%02d", 
            st->tm_year + 1900,  /* 年は +1900 が必要 */
            st->tm_mon + 1,       /* 月は +1 が必要 */
            st->tm_mday,            /* 日は 1...31 */
            st->tm_hour,            /* 時は 0...23 （0詰めの二桁で表示） */
            st->tm_min,             /* 分は 0...59 （0詰めの二桁で表示） */
            st->tm_sec);            /* 秒は 0...59 （0詰めの二桁で表示） */

    return buf;
}


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


}
#endif

