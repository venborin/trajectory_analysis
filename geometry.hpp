#include <cmath>
#include <vector>
#define R2D 57.295779513
#define D2R  0.017453293 
void help(std::string s)
{
  printf("\033[92mUsage: %s -x file -b BONDS (list) -a ANGLES (list) -d DIHEDRALS (list) -m MIXED_VARS (list)\n"
  "                  -h, --help, -?, -H: Show this help and exit.\n"
  "                  -v, --version:      print the program version and exit\n"
  "                  -s, --sep:          data separator (for writing the results)\n"
  "                  -x, --xyz:          filename.xyz          - multiple xyz file.\n"
  "                  -b, --bond:         i-j,j-k,k-l,...       - extract distances.\n"
  "                  -a, --angle:        i-j-k,k-l-m,...       - extract angles.\n"
  "                  -d, --dihedral:     i-j-k-l,j-k-l-m,...   - extract dihedrals.\n"
  "                  -o, --out-of-plane: i-j-k-l,j-k-l-m,...   - extract out-of-plane angles.\n"
  "                  -m: --mix:          i-j,k-l-m,n-o-p-q,... - extract corresponding distances, angles, dihedrals, and out-of-plane angles all together (option for lazy fucks)\033[0m\n",s.c_str());

}
void version()
{
 printf("\033[92m MD trajectory analyzer, version 1.0\n"
 " Written by Venia during the first lockdown\n"
 " Giv'at Ram, Jerusalem, Israel\n"
 " Version history:\n"
 " 1.0.0 - curent, June 22, 2021, -m mixed option for lazy assholes is added\n"
 " 0.5.1 - December 2, 2020, minor bugs are fixed\n"
 " 0.5.0 - September 17, 2020, option for calculating out-of-plane is added\n"
 " 0.0.1 - April 5, 2020, the alpha version\n\033[0m");
exit(0);
}
void split(const std::string& subject, std::vector<std::string>& container, const char *sep )
{
  container.clear();
  size_t len = subject.length() + 1;
  char* s = new char[ len ];
  memset(s, 0, len*sizeof(char));
  memcpy(s, subject.c_str(), (len - 1)*sizeof(char));
  for (char *p = strtok(s, sep); p != NULL; p = strtok(NULL, sep))
  {
    container.push_back( p );
  }
  delete[] s;
}

inline bool exists (const std::string& name) 
{
  if (FILE *file = fopen(name.c_str(), "r")) 
  {
    fclose(file);
    return true;
  } 
  else 
  {
    printf("File \"%s\" doesn't exist!\n",name.c_str());
    return false;
  }
}

void FixD(std::vector<std::vector<double> > &d)
{
for (int i = 1; i<d.size(); i++)
 {
 for (int j = 0; j < d[0].size(); j++)
  {
  if (d[i][j] - d[i-1][j] > 180) { d[i][j]-=360; }
  else if (d[i][j] - d[i-1][j] < -180) { d[i][j]+=360; }
  else d[i][j]=d[i][j];
  }
 }
}



class Geometry
{

public:
std::vector<std::vector<double>> coord;

double dist(int i, int j)
{
  return sqrt((coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0]) + (coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1]) + (coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]));
}

double unit(int c, int i, int j)
{
  return -(coord[i][c]-coord[j][c])/dist(i,j);
}

double angle(int i, int j, int k)
{
  return acos(unit(0,j,i) * unit(0,j,k) + unit(1,j,i) * unit(1,j,k) + unit(2,j,i) * unit(2,j,k))*R2D;
}

double torsion(int i, int j, int k, int l)
{
  double eabc_x = (unit(1,j,i)*unit(2,j,k) - unit(2,j,i)*unit(1,j,k));
  double eabc_y = (unit(2,j,i)*unit(0,j,k) - unit(0,j,i)*unit(2,j,k));
  double eabc_z = (unit(0,j,i)*unit(1,j,k) - unit(1,j,i)*unit(0,j,k));

  double ebcd_x = (unit(1,k,j)*unit(2,k,l) - unit(2,k,j)*unit(1,k,l));
  double ebcd_y = (unit(2,k,j)*unit(0,k,l) - unit(0,k,j)*unit(2,k,l));
  double ebcd_z = (unit(0,k,j)*unit(1,k,l) - unit(1,k,j)*unit(0,k,l));

  double exx = eabc_x * ebcd_x;
  double eyy = eabc_y * ebcd_y;
  double ezz = eabc_z * ebcd_z;

  double tau = (exx + eyy + ezz)/(sin(angle(i,j,k)*D2R) * sin(angle(j,k,l)*D2R));

  if(tau < -1.0) tau = acos(-1.0);
  else if(tau > 1.0) tau = acos(1.0);
  else tau = acos(tau);

  double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
  double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
  double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
  double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
  cross_x /= norm;
  cross_y /= norm;
  cross_z /= norm;
  double sign = 1.0;
  double dot = cross_x*unit(0,j,k)+cross_y*unit(1,j,k)+cross_z*unit(2,j,k);
  if(dot < 0.0) sign = -1.0;

  return tau*sign*R2D;
} 

double OOP(int i, int j, int k, int l)
{
  double ebcd_x = (unit(1,k,j) * unit(2,k,l) - unit(2,k,j) * unit(1,k,l));
  double ebcd_y = (unit(2,k,j) * unit(0,k,l) - unit(0,k,j) * unit(2,k,l));
  double ebcd_z = (unit(0,k,j) * unit(1,k,l) - unit(1,k,j) * unit(0,k,l));
  
  double exx = ebcd_x * unit(0,k,i);
  double eyy = ebcd_y * unit(1,k,i);
  double ezz = ebcd_z * unit(2,k,i);

  double theta = (exx + eyy + ezz)/sin(angle(j,k,l)*D2R);

  if(theta < -1.0) theta = asin(-1.0);
  else if(theta > 1.0) theta = asin(1.0);
  else theta = asin(theta);

  return theta*R2D;
}

void push(std::vector<double> &a)
{
  coord.push_back(a);
}

void clear()
{
  coord.clear();
}


};
