#include <string>
#include <cstdio>
#include <cstring>
#include "geometry.hpp"
#include <iostream>

int main(int argc, char const ** argv)
{
Geometry G;
FILE* xyz;
FILE* dihedrals;
FILE* angles;
FILE* distances;
FILE* out_of_plane;
char* sep="";
int nat;
char c[2],line[128];
float x,y,z;
size_t lineCount {0};
int frameCount=0;
bool n=false;
std::vector<double> g;
std::vector<double> tor;

std::vector<std::string> tmp;
std::vector<std::string> bnd;
std::vector<std::string> ang;
std::vector<std::string> dhd;
std::vector<std::string> oop;
std::vector<std::string> vars;
std::vector<std::string> var;
std::vector<std::vector<double> >dihds;


if (argc==1){ printf("\033[91mNo argument provided.\033[0m \n"); help(argv[0]); return 1;}

else 
{
 for(int i=1; i<argc; ++i)
 {

  std::string arg = argv[i];
  if (arg == "-h" || arg == "--help" || arg == "-?" || arg == "-H") 
  { 
  help(argv[0]);
  exit(0); 
  }
  else if (arg == "-v" || arg == "--version") {version();}
  else if (arg == "-x" || arg == "--xyz") 
  {
   if (argv[i+1] == NULL)
   {
    printf("\033[91mNo filename provided.\033[0m \n");
    help(argv[0]);
    exit(1);
   } 
   if (exists(argv[i+1])) 
   {
    xyz = fopen(argv[++i],"r");
   } 
   else 
   { 
    help(argv[0]); 
    exit(1);
   } 
  }
  else if (arg == "-b" || arg == "--bond")		{ if (argv[i+1]) { split(argv[++i],bnd,","); } else {help(argv[0]); exit(1);} } 
  else if (arg == "-a" || arg == "--angle")		{ if (argv[i+1]) { split(argv[++i],ang,","); } else {help(argv[0]); exit(1);} }
  else if (arg == "-d" || arg == "--dihedral")		{ if (argv[i+1]) { split(argv[++i],dhd,","); } else {help(argv[0]); exit(1);} }
  else if (arg == "-o" || arg == "--out-of-plane")	{ if (argv[i+1]) { split(argv[++i],oop,","); } else {help(argv[0]); exit(1);} }
  else if (arg == "-m" || arg == "--mix") 		{ if (argv[i+1]) { split(argv[++i],vars,","); } else {help(argv[0]); exit(1);} }
  else if (arg == "-s" || arg == "--sep") 		{ if (argv[i+1]) { sep=(char*)argv[++i]; } else { help(argv[0]); exit(1);} }
  else if (arg == "-e" || arg == "--enum") 		{ if (argv[i+1]) { n=true; } }
  else { printf("\033[91m Argument %s not recognized!\033[0m\n",argv[i]);help(argv[0]); exit(1); }
 }

}
if (vars.size() > 0)
{
bnd.clear();
ang.clear();
dhd.clear();
oop.clear();

for(int i=0;i < vars.size();i++)
{
 split(vars[i],var,"-");
 if (var.size() == 2)
 {
  bnd.push_back(vars[i]);
 }
 else if (var.size() == 3)
 {
  ang.push_back(vars[i]);
 }
 else if (var.size() == 4)
 {
  dhd.push_back(vars[i]);
  oop.push_back(vars[i]);
 }
 else
 {
 printf("\033[91mNo geometrical parameter is defined by %lu atoms! Exit...\033[0m\n",var.size());
 exit(1);
 }

}
}

if (bnd.size()>0 ) 
{ 
  distances = fopen("distances","w+");  
  for (int i = 0; i < bnd.size(); i++) 
  {
    fprintf(distances,"%12s%s",bnd[i].c_str(),sep);
  }
  fprintf(distances,"\n");
}
if (ang.size()>0 ) 
{ 
  angles = fopen("angles","w+"); 
  for (int i = 0; i < ang.size(); i++) 
  {
    fprintf(angles,"%12s%s",ang[i].c_str(),sep);
  }
  fprintf(angles,"\n");

}
if (dhd.size()>0 ) 
{ 
  dihedrals = fopen("dihedrals","w+"); 
  for (int i = 0; i < dhd.size(); i++)
  { 
    fprintf(dihedrals,"%12s%s",dhd[i].c_str(),sep);
  }
  fprintf(dihedrals,"\n");
}
if (oop.size() >0)
{
  out_of_plane = fopen("out_of_plane","w+");
  for (int i = 0; i < oop.size(); i++)
  {
    fprintf(out_of_plane,"%12s%s",oop[i].c_str(),sep);
  }
  fprintf(out_of_plane,"\n");
}
while (fgets(line,128,xyz)) 
{
 if (lineCount ==0 ) { sscanf(line,"%d",&nat); }
 if (lineCount%(nat+2) == 0 || (lineCount-1)%(nat+2) == 0) { ++lineCount; continue;}
 else { sscanf(line,"%s%f%f%f",c,&x,&y,&z); g.push_back(x); g.push_back(y);g.push_back(z); G.push(g); g.clear();++lineCount;}
 if (lineCount%(nat+2) == 0) 
 { 
 if (bnd.size()>0)
  {
  if (n==true) { fprintf(distances,"%5d",frameCount);}
  for (int i = 0; i < bnd.size(); i++)
  {
   split(bnd[i],tmp,"-");
   fprintf(distances,"%12.6f%s",G.dist(stoi(tmp[0])-1,stoi(tmp[1])-1),sep);
   tmp.clear();
  }
  fprintf(distances,"\n");
  }
 if (ang.size() >0)
  {
  for (int i = 0; i < ang.size(); i++)
  {
   split(ang[i],tmp,"-");
   fprintf(angles,"%12.6f%s",G.angle(stoi(tmp[0])-1,stoi(tmp[1])-1,stoi(tmp[2])-1),sep);
   tmp.clear();
  }
  fprintf(angles,"\n");
  }
 if (dhd.size() >0)
  {
  for (int i = 0; i < dhd.size(); i++)
  {
   split(dhd[i],tmp,"-");
   tor.push_back(G.torsion(stoi(tmp[0])-1,stoi(tmp[1])-1,stoi(tmp[2])-1,stoi(tmp[3])-1));
   tmp.clear();
  }
  dihds.push_back(tor); tor.clear();
  }
 if (oop.size() >0)
  {
  for (int i = 0; i < oop.size(); i++)
  {
   split(oop[i],tmp,"-");
   fprintf(out_of_plane,"%12.6f%s",G.OOP(stoi(tmp[0])-1,stoi(tmp[1])-1,stoi(tmp[2])-1,stoi(tmp[3])-1),sep);
   tmp.clear();
  }
  fprintf(out_of_plane,"\n"); 
  }
frameCount++;
G.clear();
}
}
if (ang.size() >0) { fclose(angles); }
if (bnd.size() >0) { fclose(distances); }
if (oop.size() >0) { fclose(out_of_plane); }
if (dhd.size() >0) 
{
  FixD(dihds);
  for (int i = 0; i < dihds.size(); i++) {
    for (int j = 0; j < dihds[0].size(); j++) {
      fprintf(dihedrals,"%12.6f%s",dihds[i][j],sep);
    }
    fprintf(dihedrals,"\n");
  }
  fclose(dihedrals);
}


fclose(xyz);

}
