#pragma once
#include <cstdint>
#include <cmath>
#include <limits>

#include "d3dx9math.h"

namespace IKD
{

    namespace aux{
        template<class T>
        struct vector{
        public:
            union{
                struct{
                    T x, y, z;
                };
                T coords[3];
            };

            vector() = default;

            vector(const T *ptr){
                x = ptr[0];
                y = ptr[1];
                z = ptr[2];
            }

            vector(const vector &other) : x(other.x), y(other.y), z(other.z){}

            vector(T x, T y, T z) : x(x), y(y), z(z){}

            operator T*(){
                return &coords[0];
            }

            operator const T*() const{
                return &coords[0];
            }

            vector &operator =(const vector &other){
                for(int i = 0; i < 3; ++i){
                    coords[i] = other.coords[i];
                }
                return *this;
            }

            vector &operator +=(const vector &other){
                for(int i = 0; i < 3; ++i){
                    coords[i] += other.coords[i];
                }
                return *this;
            }

            vector &operator -=(const vector &other){
                for(int i = 0; i < 3; ++i){
                    coords[i] -= other.coords[i];
                }
                return *this;
            }

            vector &operator *=(T a){
                for(int i = 0; i < 3; ++i){
                    coords[i] *= a;
                }
                return *this;
            }

            vector &operator /=(T a){
                for(int i = 0; i < 3; ++i){
                    coords[i] /= a;
                }
                return *this;
            }

            vector operator +() const{
                return *this;
            }

            vector operator -() const{
                vector other(*this);
                for(int i = 0; i < 3; ++i){
                    other.coords[i] = -other.coords[i];
                }
                return other;
            }

            vector operator +(const vector &other) const{
                vector ret(*this);
                ret += other;
                return ret;
            }

            vector operator -(const vector &other) const{
                vector ret(*this);
                ret -= other;
                return ret;
            }

            vector operator *(T a) const{
                vector ret(*this);
                ret *= a;
                return ret;
            }

            vector operator /(T a) const{
                vector ret(*this);
                ret /= a;
                return ret;
            }

            bool operator ==(const vector &other) const{
                for(int i = 0; i < 3; ++i){
                    if(coords[i] != other.coords[i]){
                        return false;
                    }
                }
                return true;
            }

            bool operator !=(const vector &other) const{
                return !(*this == other);
            }

            static T LengthSq(const vector *ptr){
                T r2 = T(0);
                for(int i = 0; i < 3; ++i){
                    r2 += ptr->coords[i] * ptr->coords[i];
                }
                return r2;
            }

            static T Vec3Dot(const vector *a, const vector *b){
                T d = T(0);
                for(int i = 0; i < 3; ++i){
                    d += a->coords[i] * b->coords[i];
                }
                return d;
            }

            static vector *Vec3Normalize(vector *out, const vector *in){
                *out = *in / std::sqrt(Vec3Dot(in, in));
                return out;
            }
        };
    }

    using vector = aux::vector<double>;

///////////////////////////////////////////////////
// パーティクル衝突判定・時刻・位置算出関数
//   rA          : パーティクルAの半径
//   rB          : パーティクルBの半径
//   pre_pos_A   : パーティクルAの前の位置
//   pos_A       : パーティクルAの次の到達位置
//   pre_pos_B   : パーティクルBの前位置
//   pos_B       : パーティクルBの次の到達位置
//   pout_t      : 衝突時間を格納するdouble型へのポインタ
//   pout_colli_A : パーティクルAの衝突位置を格納するD3DXVECTOR型へのポインタ
//   pout_colli_B : パーティクルAの衝突位置を格納するD3DXVECTOR型へのポインタ

bool CalcParticleCollision(
   double rA, double rB, 
   const vector &pPre_pos_A, const vector &pPos_A,
   const vector &pPre_pos_B, const vector &pPos_B,
   double &pOut_t,
   vector &pOut_colli_A,
   vector &pOut_colli_B
)
{
   // 前位置及び到達位置におけるパーティクル間のベクトルを算出
   vector C0 = pPre_pos_B - pPre_pos_A;
   vector C1 = pPos_B - pPos_A;
   vector D = C1 - C0;

   // 衝突判定用の2次関数係数の算出
   double P = vector::LengthSq( &D ); if(P==0) return false; // 同じ方向に移動
   double Q = vector::Vec3Dot( &C0, &D );
   double R = vector::LengthSq( &C0 );

   // パーティクル距離
   double r = rA + rB;

   // 衝突判定式
   double Judge = Q*Q - P*(R-r*r);
   if( Judge < 0 ){
      // 衝突していない
      return false;
   }

   // 衝突時間の算出
   double t_plus = (-Q + std::sqrt(Judge))/P;
   double t_minus = (-Q - std::sqrt(Judge))/P;

   // 衝突時間が0未満1より大きい場合、衝突しない
   if( (t_plus < 0 || t_plus > 1) && (t_minus < 0 || t_minus > 1)) return false;

   // 衝突時間の決定（t_minus側が常に最初の衝突）
   pOut_t = t_minus;

   // 衝突位置の決定
   pOut_colli_A = pPre_pos_A + (pPos_A - pPre_pos_A) * t_minus;
   pOut_colli_B = pPre_pos_B + (pPos_B - pPre_pos_B) * t_minus;

   return true; // 衝突報告
}




///////////////////////////////////////////////////
// パーティクル衝突後速度位置算出関数
//   pColliPos_A : 衝突中のパーティクルAの中心位置
//   pVelo_A     : 衝突の瞬間のパーティクルAの速度
//   pColliPos_B : 衝突中のパーティクルBの中心位置
//   pVelo_B     : 衝突の瞬間のパーティクルBの速度
//   weight_A    : パーティクルAの質量
//   weight_B    : パーティクルBの質量
//   res_A       : パーティクルAの反発率
//   res_B       : パーティクルBの反発率
//   time        : 反射後の移動可能時間
//   pOut_pos_A  : パーティクルAの反射後位置
//   pOut_velo_A : パーティクルAの反射後速度ベクトル
//   pOut_pos_B  : パーティクルBの反射後位置
//   pOut_velo_B : パーティクルBの反射後速度ベクトル
bool CalcParticleColliAfterPos(
   const vector &pColliPos_A, const vector &pVelo_A,
    const vector &pColliPos_B, const vector &pVelo_B,
   double weight_A, double weight_B,
   double res_A, double res_B,
   double time,
   vector &pOut_pos_A, vector &pOut_velo_A,
   vector &pOut_pos_B, vector &pOut_velo_B
)
{
   double TotalWeight = weight_A + weight_B; // 質量の合計
   double RefRate = (1 + res_A*res_B); // 反発率
   vector C = pColliPos_B - pColliPos_A; // 衝突軸ベクトル
   vector::Vec3Normalize(&C, &C);
   double Dot = vector::Vec3Dot( &(pVelo_A-pVelo_B), &C ); // 内積算出
   vector ConstVec = C * (RefRate*Dot/TotalWeight); // 定数ベクトル

   // 衝突後速度ベクトルの算出
   pOut_velo_A = ConstVec * -weight_B + pVelo_A;
   pOut_velo_B = ConstVec * weight_A + pVelo_B;

   // 衝突後位置の算出
   pOut_pos_A = pColliPos_A + (pOut_velo_A) * time;
   pOut_pos_B = pColliPos_B + (pOut_velo_B) * time;

   return true;
}


#define COLLISION_EPSIRON 0.00001 // 誤差

///////////////////////////////////////////////////
// 平面パーティクル衝突判定・時刻・位置算出関数
// r : パーティクルの半径
// pPre_pos : パーティクルの前の位置
// pPos : パーティクルの次の到達位置
// pNormal : 平面の法線
// pPlane_pos : 平面上の1点
// pOut_t : 衝突時間を格納するdouble型へのポインタ
// pOut_colli : パーティクルの衝突位置を格納するD3DXVECTOR型へのポインタ
// 戻り値 : 衝突(true), 非衝突(false)

bool CalcParticlePlaneCollision(
   double r,
   const vector &pPre_pos, const vector &pPos,
    const vector &pNormal, const vector &pPlane_pos,
   double &t,
   vector &pOut_colli
) 
{
   vector C0 = pPre_pos - pPlane_pos; // 平面上の一点から現在位置へのベクトル
   vector D = pPos - pPre_pos; // 現在位置から予定位置までのベクトル
   vector N; // 法線
   vector::Vec3Normalize(&N, &pNormal); // 法線を標準化

   // 平面と中心点の距離を算出
   double Dot_C0 = vector::Vec3Dot( &C0, &N );
   double dist_plane_to_point = std::fabs( Dot_C0 );

   // 進行方向と法線の関係をチェック
   double Dot = vector::Vec3Dot( &D, &N );

   // 平面と平行に移動してめり込んでいるスペシャルケース
   if( (COLLISION_EPSIRON-std::fabs(Dot) > 0.0) && (dist_plane_to_point < r) ){
      // 一生抜け出せないので最大時刻を返す
      t = (std::numeric_limits<double>::max)();
      // 衝突位置は仕方ないので今の位置を指定
      pOut_colli = pPre_pos;
      return true;
   }

   // 交差時間の算出
   t = ( r - Dot_C0 )/Dot;

   // 衝突位置の算出
   pOut_colli = pPre_pos + D * t;

   // めり込んでいたら衝突として処理終了
   if ( dist_plane_to_point < r )
      return true;

   // 壁に対して移動が逆向きなら衝突していない
   if( Dot >= 0 )
      return false;

   // 時間が0～1の間にあれば衝突
   if( (0 <= t) && (t <= 1) )
      return true;

   return false;
}


///////////////////////////////////////////////////
// 平面と球の衝突後速度算出関数
// pColliPos : 衝突中のパーティクルの中心位置
// pVelo : 衝突の瞬間のパーティクルの速度
// res : パーティクルの壁に対する反発率
// time : 反射後の移動可能時間
// pNormal : 平面の法線
// pOut_pos : パーティクルの反射後位置
// pOut_velo : パーティクルの反射後速度ベクトル

bool CalcParticlePlaneAfterPos(
    const vector &pColliPos,
    const vector &pVelo,
    double res,
    double time,
    vector &pNormal,
    vector &pOut_pos,
    vector &pOut_velo
)
{
    // 反射後速度を算出
    vector N;
    vector::Vec3Normalize(&N,&pNormal);
    pOut_velo = pVelo - N * ((1+res)*vector::Vec3Dot(&N,&pVelo));

    // 移動位置を計算
    pOut_pos = pColliPos + pOut_velo * time;

    return true;
}

#undef COLLISION_EPSIRON


}    // end namespace IKD