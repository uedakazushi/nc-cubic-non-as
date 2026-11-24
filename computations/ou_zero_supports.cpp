// Necessary support conditions for the cyclic three-block quiver.
// Bit 9*i+3*j+k denotes a nonzero coefficient c_ijk.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <vector>
using namespace std;
int root(int *p,int x){while(p[x]!=x)x=p[x];return x;}
void join(int*p,int a,int b){p[root(p,a)]=root(p,b);}
bool connected(int a){int p[6]={0,1,2,3,4,5};for(int j=0;j<3;j++)for(int k=0;k<3;k++)if(a>>(3*j+k)&1)join(p,j,k+3);for(int j=1;j<6;j++)if(root(p,j)!=root(p,0))return false;return true;}
bool ordered_pair_ok(int a,int b){int p[6]={0,1,2,3,4,5},force[6]={},verts[6]={},edges[6]={},f[6]={};for(int j=0;j<3;j++)for(int k=0;k<3;k++){int x=1<<(3*j+k);if((a&b)&x)join(p,j,k+3);else if(a&x)force[k+3]=1;else if(b&x)force[j]=1;}for(int j=0;j<6;j++){int r=root(p,j);verts[r]++;f[r]|=force[j];}for(int j=0;j<3;j++)for(int k=0;k<3;k++)if((a&b)>>(3*j+k)&1)edges[root(p,j)]++;for(int j=0;j<6;j++)if(verts[j]&&!f[j]&&edges[j]<verts[j])return false;return true;}
int slice(uint32_t x,int d,int a){int m=0;for(int i=0;i<3;i++)for(int j=0;j<3;j++){int t[3];t[d]=a;t[(d+1)%3]=i;t[(d+2)%3]=j;if(x>>(9*t[0]+3*t[1]+t[2])&1)m|=1<<(3*i+j);}return m;}
int main(){
 vector<int> v; static bool ok[512][512]={}; bool conn[512]={};
 for(int a=0;a<512;a++)if(connected(a)){v.push_back(a);conn[a]=true;}
 assert(v.size()==205);
 for(int a:v)for(int b:v)ok[a][b]=ordered_pair_ok(a,b)&&ordered_pair_ok(b,a);
 vector<uint32_t> candidates; long long ct[13]={},cnt[13]={};
 for(int a:v)for(int b:v)for(int c:v){uint32_t t=a|(b<<9)|(c<<18);int n=27-__builtin_popcount(t);int z[2][3];bool good=true;for(int d=1;d<3&&good;d++)for(int i=0;i<3;i++){z[d-1][i]=slice(t,d,i);if(!conn[z[d-1][i]]){good=false;break;}}if(!good)continue;ct[n]++;if(!ok[a][b]||!ok[a][c]||!ok[b][c])continue;for(int d=0;d<2&&good;d++)for(int i=0;i<3;i++)for(int j=i+1;j<3;j++)if(!ok[z[d][i]][z[d][j]])good=false;if(!good)continue;cnt[n]++;candidates.push_back(t);}
 vector<array<int,3>> ps;array<int,3> p={0,1,2};do{ps.push_back(p);}while(next_permutation(p.begin(),p.end()));
 vector<array<int,27>> maps;
 for(auto ax:ps)for(auto a:ps)for(auto b:ps)for(auto c:ps){array<int,27> m;array<array<int,3>,3> per={a,b,c};for(int x=0;x<27;x++){int t[]={x/9,x/3%3,x%3};m[x]=9*per[0][t[ax[0]]]+3*per[1][t[ax[1]]]+per[2][t[ax[2]]];}maps.push_back(m);}
 vector<bool> seen(1<<27,false);int counts[13]={};
 for(uint32_t x:candidates){if(seen[x])continue;int n=27-__builtin_popcount(x);counts[n]++;cout<<n<<" "<<x<<"\n";for(auto &m:maps){uint32_t y=0;for(int i=0;i<27;i++)if(x>>i&1)y|=1u<<m[i];seen[y]=true;}}
 const int exp_orbits[]={1,1,3,9,30,80,201,344,387,217,53,4,1};
 const long long exp_counts[]={1,27,351,2898,16659,68256,196110,374490,428004,242688,49680,1944,72};
 for(int n=0;n<=12;n++){assert(counts[n]==exp_orbits[n]);assert(cnt[n]==exp_counts[n]);cerr<<n<<" "<<ct[n]<<" "<<cnt[n]<<" "<<counts[n]<<"\n";}
}
