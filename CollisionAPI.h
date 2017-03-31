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
// �p�[�e�B�N���Փ˔���E�����E�ʒu�Z�o�֐�
//   rA          : �p�[�e�B�N��A�̔��a
//   rB          : �p�[�e�B�N��B�̔��a
//   pre_pos_A   : �p�[�e�B�N��A�̑O�̈ʒu
//   pos_A       : �p�[�e�B�N��A�̎��̓��B�ʒu
//   pre_pos_B   : �p�[�e�B�N��B�̑O�ʒu
//   pos_B       : �p�[�e�B�N��B�̎��̓��B�ʒu
//   pout_t      : �Փˎ��Ԃ��i�[����double�^�ւ̃|�C���^
//   pout_colli_A : �p�[�e�B�N��A�̏Փˈʒu���i�[����D3DXVECTOR�^�ւ̃|�C���^
//   pout_colli_B : �p�[�e�B�N��A�̏Փˈʒu���i�[����D3DXVECTOR�^�ւ̃|�C���^

bool CalcParticleCollision(
   double rA, double rB, 
   const vector &pPre_pos_A, const vector &pPos_A,
   const vector &pPre_pos_B, const vector &pPos_B,
   double &pOut_t,
   vector &pOut_colli_A,
   vector &pOut_colli_B
)
{
   // �O�ʒu�y�ѓ��B�ʒu�ɂ�����p�[�e�B�N���Ԃ̃x�N�g�����Z�o
   vector C0 = pPre_pos_B - pPre_pos_A;
   vector C1 = pPos_B - pPos_A;
   vector D = C1 - C0;

   // �Փ˔���p��2���֐��W���̎Z�o
   double P = vector::LengthSq( &D ); if(P==0) return false; // ���������Ɉړ�
   double Q = vector::Vec3Dot( &C0, &D );
   double R = vector::LengthSq( &C0 );

   // �p�[�e�B�N������
   double r = rA + rB;

   // �Փ˔��莮
   double Judge = Q*Q - P*(R-r*r);
   if( Judge < 0 ){
      // �Փ˂��Ă��Ȃ�
      return false;
   }

   // �Փˎ��Ԃ̎Z�o
   double t_plus = (-Q + std::sqrt(Judge))/P;
   double t_minus = (-Q - std::sqrt(Judge))/P;

   // �Փˎ��Ԃ�0����1���傫���ꍇ�A�Փ˂��Ȃ�
   if( (t_plus < 0 || t_plus > 1) && (t_minus < 0 || t_minus > 1)) return false;

   // �Փˎ��Ԃ̌���it_minus������ɍŏ��̏Փˁj
   pOut_t = t_minus;

   // �Փˈʒu�̌���
   pOut_colli_A = pPre_pos_A + (pPos_A - pPre_pos_A) * t_minus;
   pOut_colli_B = pPre_pos_B + (pPos_B - pPre_pos_B) * t_minus;

   return true; // �Փ˕�
}




///////////////////////////////////////////////////
// �p�[�e�B�N���Փˌ㑬�x�ʒu�Z�o�֐�
//   pColliPos_A : �Փ˒��̃p�[�e�B�N��A�̒��S�ʒu
//   pVelo_A     : �Փ˂̏u�Ԃ̃p�[�e�B�N��A�̑��x
//   pColliPos_B : �Փ˒��̃p�[�e�B�N��B�̒��S�ʒu
//   pVelo_B     : �Փ˂̏u�Ԃ̃p�[�e�B�N��B�̑��x
//   weight_A    : �p�[�e�B�N��A�̎���
//   weight_B    : �p�[�e�B�N��B�̎���
//   res_A       : �p�[�e�B�N��A�̔�����
//   res_B       : �p�[�e�B�N��B�̔�����
//   time        : ���ˌ�̈ړ��\����
//   pOut_pos_A  : �p�[�e�B�N��A�̔��ˌ�ʒu
//   pOut_velo_A : �p�[�e�B�N��A�̔��ˌ㑬�x�x�N�g��
//   pOut_pos_B  : �p�[�e�B�N��B�̔��ˌ�ʒu
//   pOut_velo_B : �p�[�e�B�N��B�̔��ˌ㑬�x�x�N�g��
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
   double TotalWeight = weight_A + weight_B; // ���ʂ̍��v
   double RefRate = (1 + res_A*res_B); // ������
   vector C = pColliPos_B - pColliPos_A; // �Փˎ��x�N�g��
   vector::Vec3Normalize(&C, &C);
   double Dot = vector::Vec3Dot( &(pVelo_A-pVelo_B), &C ); // ���ώZ�o
   vector ConstVec = C * (RefRate*Dot/TotalWeight); // �萔�x�N�g��

   // �Փˌ㑬�x�x�N�g���̎Z�o
   pOut_velo_A = ConstVec * -weight_B + pVelo_A;
   pOut_velo_B = ConstVec * weight_A + pVelo_B;

   // �Փˌ�ʒu�̎Z�o
   pOut_pos_A = pColliPos_A + (pOut_velo_A) * time;
   pOut_pos_B = pColliPos_B + (pOut_velo_B) * time;

   return true;
}


#define COLLISION_EPSIRON 0.00001 // �덷

///////////////////////////////////////////////////
// ���ʃp�[�e�B�N���Փ˔���E�����E�ʒu�Z�o�֐�
// r : �p�[�e�B�N���̔��a
// pPre_pos : �p�[�e�B�N���̑O�̈ʒu
// pPos : �p�[�e�B�N���̎��̓��B�ʒu
// pNormal : ���ʂ̖@��
// pPlane_pos : ���ʏ��1�_
// pOut_t : �Փˎ��Ԃ��i�[����double�^�ւ̃|�C���^
// pOut_colli : �p�[�e�B�N���̏Փˈʒu���i�[����D3DXVECTOR�^�ւ̃|�C���^
// �߂�l : �Փ�(true), ��Փ�(false)

bool CalcParticlePlaneCollision(
   double r,
   const vector &pPre_pos, const vector &pPos,
    const vector &pNormal, const vector &pPlane_pos,
   double &t,
   vector &pOut_colli
) 
{
   vector C0 = pPre_pos - pPlane_pos; // ���ʏ�̈�_���猻�݈ʒu�ւ̃x�N�g��
   vector D = pPos - pPre_pos; // ���݈ʒu����\��ʒu�܂ł̃x�N�g��
   vector N; // �@��
   vector::Vec3Normalize(&N, &pNormal); // �@����W����

   // ���ʂƒ��S�_�̋������Z�o
   double Dot_C0 = vector::Vec3Dot( &C0, &N );
   double dist_plane_to_point = std::fabs( Dot_C0 );

   // �i�s�����Ɩ@���̊֌W���`�F�b�N
   double Dot = vector::Vec3Dot( &D, &N );

   // ���ʂƕ��s�Ɉړ����Ă߂荞��ł���X�y�V�����P�[�X
   if( (COLLISION_EPSIRON-std::fabs(Dot) > 0.0) && (dist_plane_to_point < r) ){
      // �ꐶ�����o���Ȃ��̂ōő厞����Ԃ�
      t = (std::numeric_limits<double>::max)();
      // �Փˈʒu�͎d���Ȃ��̂ō��̈ʒu���w��
      pOut_colli = pPre_pos;
      return true;
   }

   // �������Ԃ̎Z�o
   t = ( r - Dot_C0 )/Dot;

   // �Փˈʒu�̎Z�o
   pOut_colli = pPre_pos + D * t;

   // �߂荞��ł�����Փ˂Ƃ��ď����I��
   if ( dist_plane_to_point < r )
      return true;

   // �ǂɑ΂��Ĉړ����t�����Ȃ�Փ˂��Ă��Ȃ�
   if( Dot >= 0 )
      return false;

   // ���Ԃ�0�`1�̊Ԃɂ���ΏՓ�
   if( (0 <= t) && (t <= 1) )
      return true;

   return false;
}


///////////////////////////////////////////////////
// ���ʂƋ��̏Փˌ㑬�x�Z�o�֐�
// pColliPos : �Փ˒��̃p�[�e�B�N���̒��S�ʒu
// pVelo : �Փ˂̏u�Ԃ̃p�[�e�B�N���̑��x
// res : �p�[�e�B�N���̕ǂɑ΂��锽����
// time : ���ˌ�̈ړ��\����
// pNormal : ���ʂ̖@��
// pOut_pos : �p�[�e�B�N���̔��ˌ�ʒu
// pOut_velo : �p�[�e�B�N���̔��ˌ㑬�x�x�N�g��

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
    // ���ˌ㑬�x���Z�o
    vector N;
    vector::Vec3Normalize(&N,&pNormal);
    pOut_velo = pVelo - N * ((1+res)*vector::Vec3Dot(&N,&pVelo));

    // �ړ��ʒu���v�Z
    pOut_pos = pColliPos + pOut_velo * time;

    return true;
}

#undef COLLISION_EPSIRON


}    // end namespace IKD