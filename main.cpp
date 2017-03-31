#include "CollisionAPI.h"    // �Փ˔���API
#include "ColTrees.h"    // ���`4���؃w�b�_�[
#include <cstdint>

#include <DxLib.h>
#include <windows.h>

using namespace std;
using namespace IKD;


///////////////////////////////////////////////
// �O���[�o���ϐ��Q�i���l��ς��ėV�ׂ܂��j
////////////////
double g_Circle_Ref = 0.90;    // �����m�̔����W��
double Wall_Ref = 0.90;        // ���ƕǂ̔����W��
const int g_CircleNum = 1500;    // ��������~�̐�
int g_PartitionLebel = 7;        // ��ԕ������x��
double g_Gravity = 9.8;            // �d��


// �~�\����
struct CIRCLE
{
    std::uint32_t ID;
    double x, y;                // �~�̈ʒu
    double Pre_x, Pre_y;        // 1�O�̉~�̈ʒu
    double vx, vy;            // ���x�x�N�g��
    double ax, ay;            // �����x�x�N�g��
    double r;                // ���a
    double w;                // ����
    double scale;            // �X�P�[��

    CIRCLE()
    {
        x = y = vx = vy = ax = ay = Pre_x = Pre_y = 0.0f;
        r = 1.0f;
        w = 1.0f;
        scale = 1.0f;
    }
};


// �ǂƋ��̔��˃x�N�g�����v�Z
void GetRefrectVelo(IKD::vector *pOut, IKD::vector &N, IKD::vector &V, double e)
{
    IKD::vector::Vec3Normalize(&N,&N);
    *pOut = V - N*((1+e)*IKD::vector::Vec3Dot(&N,&V));
}


// �ǂƂ̔��ˌ�̈ʒu���Z�o
void GetRelectedPos( double Res_time, CIRCLE &circle, IKD::vector &RefV )
{
    // �Փˈʒu
    // 0.99�͕ǂɂ߂荞�܂Ȃ����߂̕␳
    circle.x = circle.Pre_x + circle.vx * (1-Res_time)*0.99f;
    circle.y = circle.Pre_y + circle.vy * (1-Res_time)*0.99f;
    // ���˃x�N�g��
    circle.vx = RefV.x;
    circle.vy = RefV.y;
    // �ʒu��␳
    circle.x += circle.vx * Res_time;
    circle.y += circle.vy * Res_time;
}


// ���̉~�̈ʒu���Z�o
void GetNextCirclePos( CIRCLE &circle )
{
    IKD::vector RefV;    // ���ˌ�̑��x�x�N�g��
    double Res_time = 0.0f;    // �Փˌ�̈ړ��\����

    // �d�͂��|���ė��Ƃ�
    circle.vy += g_Gravity/60;    // 1�t���[����9.8/60(m/s)����

    // ���̑��x�ňʒu���X�V
    circle.Pre_x = circle.x;    // �O�̈ʒu��ۑ�
    circle.Pre_y = circle.y;
    circle.x += circle.vx;        // �ʒu�X�V
    circle.y += circle.vy;

    // �ǂƂ̏Փ˂��`�F�b�N
    // ����
    if(circle.x<0){
        // ���ˌ�̑��x�x�N�g�����擾
        GetRefrectVelo( &RefV, IKD::vector(1,0,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        // �c�莞�ԎZ�o
        Res_time = circle.x / circle.vx;
        // ���ˌ�̈ʒu���Z�o
        GetRelectedPos(    Res_time, circle, RefV );
    }
    // �E��
    else if(circle.x>640){
        GetRefrectVelo( &RefV, IKD::vector(-1,0,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        Res_time = (circle.x-640) / circle.vx;
        GetRelectedPos(    Res_time, circle, RefV );
    }
    // ����
    else if(circle.y>480){
        GetRefrectVelo( &RefV, IKD::vector(0,-1,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        Res_time = (circle.y-480) / circle.vy;
        GetRelectedPos(    Res_time, circle, RefV );
    }
}


// 2�~�Փˏ���
void CircleColProc( CIRCLE *c1, CIRCLE *c2 )
{
    double t=0;
    IKD::vector C1ColPos, C2ColPos, C1Velo, C2Velo;

    // �Փ˂��Ă���2�~�̏Փˈʒu�����o
    if(!CalcParticleCollision(c1->r, c2->r,
        IKD::vector(c1->Pre_x, c1->Pre_y,0),
        IKD::vector(c1->x, c1->y, 0),
        IKD::vector(c2->Pre_x, c2->Pre_y, 0),
        IKD::vector(c2->x, c2->y, 0),
        t,
        C1ColPos,
        C2ColPos))
        return;    // �Փ˂��Ă��Ȃ��悤�ł�

    // �Փˈʒu��O�ʒu�Ƃ��ĕۑ�
    c1->x = C1ColPos.x;
    c1->y = C1ColPos.y; 
    c2->x = C2ColPos.x;
    c2->y = C2ColPos.y; 
    c1->Pre_x = C1ColPos.x;
    c1->Pre_y = C1ColPos.y; 
    c2->Pre_x = C2ColPos.x;
    c2->Pre_y = C2ColPos.y; 

    // �Փˌ�̑��x���Z�o
    if(!CalcParticleColliAfterPos(
        C1ColPos, IKD::vector(c1->vx, c1->vy, 0),
        C2ColPos, IKD::vector(c2->vx, c2->vy, 0),
        c1->w, c2->w,
        g_Circle_Ref, g_Circle_Ref,        // ���̔����W��
        t,
        C1ColPos, C1Velo,
        C2ColPos, C2Velo))
        return; // �������s�����悤�ł�

    // �Փˌ�ʒu�Ɉړ�
    c1->vx = C1Velo.x;
    c1->vy = C1Velo.y;
    c2->vx = C2Velo.x;
    c2->vy = C2Velo.y;
    c1->x += c1->vx;
    c1->y += c1->vy;
    c2->x += c2->vx;
    c2->y += c2->vy;
}

// ���C���֐�
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPTSTR lpCmdLine, int nCmdShow){
    SetOutApplicationLogValidFlag(FALSE);
    ChangeWindowMode(TRUE);
    SetGraphMode(640, 480, 32);
    SetMainWindowText("2d collision");
    if(DxLib_Init() != 0){
        return -1;
    }

    int main_graphic_handle = MakeScreen(640, 480);

    srand(0x09485524);    // �����ď�����

    //WaitKey();

    int ff = sizeof(cell<CIRCLE>);

    /////////////////////////////
    // �~�I�u�W�F�N�g�̏�����
    CIRCLE CAry[ g_CircleNum ];
    tree_object<CIRCLE> *spOFTAry[g_CircleNum];    // �~�I�u�W�F�N�g����OFT�I�u�W�F�N�g

    //////////////////////////////////////////
    // �~�I�u�W�F�N�g�̏����ʒu�E���x�̐ݒ�
    //  �V�ԃ|�C���g�ł�(^-^)
    //////////
    int cn;
    for(cn=0;cn<g_CircleNum;cn++)    // g_CircleNum�����~�𐶐�
    {
        CAry[cn].ID = cn;
        CAry[cn].r = 2 + 3*(double)rand()/RAND_MAX;    // �~�̔��a(2�`5�܂Ń����_��)
        CAry[cn].x = (double)cn/g_CircleNum*120+30*(double)rand()/RAND_MAX;        // �ג��������ʒu
        CAry[cn].y = -400 + 700*(double)cn/g_CircleNum;    // ���\�����Ƃ��납�痎�Ƃ��܂�
        CAry[cn].vx = 0.5;        // �����i�K���j
        CAry[cn].vy = 0;
        CAry[cn].scale = CAry[cn].r/32.0f;                // �摜��64�~64�Ȃ̂ŃX�P�[���l�͂����Ȃ��ł��I
        CAry[cn].w = CAry[cn].r*CAry[cn].r*CAry[cn].r;    // ���ʂ͔��a��3��ɔ��Ƃ��܂�
        // OFT�ɓo�^
        tree_object<CIRCLE> *p = new tree_object<CIRCLE>( cn );
        p->object = &CAry[cn];    // �o�^
        spOFTAry[cn] = p;
    }


    ///////////////////////////////////
    // ���`4���؃}�l�[�W��
    //  ��Ԕ͈͂�X=-60�`720; Y=-1200�`520�ɐݒ�
    //  �~����яo���Ȃ��͈͂��w�肷��Ηǂ��̂�
    //  �A�o�E�g�ł�
    liner_for_tree_manager<CIRCLE> LTree;
    if(!LTree.init(g_PartitionLebel, -60.0f, -1200.0f, 720.0f, 520.0f))
    {
        return -1;
    }

    // ���[�v���ꎞ�ϐ�
    std::uint32_t ColNum;                // �Փ˔����
    collision_list<CIRCLE>* ColVect;    // �ՓˑΏۃ��X�g

    while(true){
        Sleep(16);
        if(ProcessMessage() != 0){
            break;
        }

        SetDrawScreen(main_graphic_handle);
        SetDrawScreen(DX_SCREEN_BACK);
        DrawBox(0, 0, 640, 480, GetColor(0xFF, 0xFF, 0xFF), TRUE);

        // �~�̈ʒu���X�V���ăc���[�ɍēo�^
        static bool change = true;
        for( cn = 0; cn < g_CircleNum; cn++ )
        {
            int i = cn;
            // ���e�X�g�̂��߂Ɉ�E�ēo�^�̏��Ԃ�ς��Ă��܂�
            if ( change == true ) {
                i = g_CircleNum - cn - 1;
            }
            CIRCLE *pTmp = spOFTAry[i]->object;
            GetNextCirclePos( *pTmp );    // ���̈ړ��ʒu��������
            spOFTAry[i]->remove();        // ��x���X�g����O���
            // �ēo�^
            LTree.register_object( pTmp->x-pTmp->r, pTmp->y-pTmp->r, pTmp->x+pTmp->r, pTmp->y+pTmp->r, spOFTAry[i] );
        }
        change = !change;

        // �ՓˑΉ����X�g���擾
        ColNum = LTree.get_all_collision_list( ColVect );

        // �Փ˔���
        std::uint32_t c;
        ColNum/=2;    // 2�Ŋ���̂̓y�A�ɂȂ��Ă���̂�
        CIRCLE** pRoot = ColVect->root();
        for(c=0; c<ColNum; c++){
            double r2 = (pRoot[c*2]->r+pRoot[c*2+1]->r)*(pRoot[c*2]->r+pRoot[c*2+1]->r);
            double x = (pRoot[c*2]->x-pRoot[c*2+1]->x);
            double y = (pRoot[c*2]->y-pRoot[c*2+1]->y);
            if(r2 >= x*x + y*y )
            {
                // 2�~�Փˏ���������
                CircleColProc( pRoot[c*2], pRoot[c*2+1] );
            }
        }

        // �`��
        for(int i = 0; i < g_CircleNum; ++i){
            DrawCircle(
                static_cast<int>(CAry[i].x),
                static_cast<int>(CAry[i].y),
                static_cast<int>(CAry[i].scale * 32),
                GetColor(0, 0, 0),
                TRUE
            );
        }
        
        ScreenFlip();
    }

    for(int i = 0; i < g_CircleNum; ++i){
        delete spOFTAry[i];
    }

    DxLib_End();
   return 0;
}