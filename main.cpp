#include "CollisionAPI.h"    // 衝突判定API
#include "ColTrees.h"    // 線形4分木ヘッダー
#include <cstdint>

#include <DxLib.h>
#include <windows.h>

using namespace std;
using namespace IKD;


///////////////////////////////////////////////
// グローバル変数群（数値を変えて遊べます）
////////////////
double g_Circle_Ref = 0.90;    // 球同士の反発係数
double Wall_Ref = 0.90;        // 球と壁の反発係数
const int g_CircleNum = 1500;    // 生成する円の数
int g_PartitionLebel = 7;        // 空間分割レベル
double g_Gravity = 9.8;            // 重力


// 円構造体
struct CIRCLE
{
    std::uint32_t ID;
    double x, y;                // 円の位置
    double Pre_x, Pre_y;        // 1つ前の円の位置
    double vx, vy;            // 速度ベクトル
    double ax, ay;            // 加速度ベクトル
    double r;                // 半径
    double w;                // 質量
    double scale;            // スケール

    CIRCLE()
    {
        x = y = vx = vy = ax = ay = Pre_x = Pre_y = 0.0f;
        r = 1.0f;
        w = 1.0f;
        scale = 1.0f;
    }
};


// 壁と球の反射ベクトルを計算
void GetRefrectVelo(IKD::vector *pOut, IKD::vector &N, IKD::vector &V, double e)
{
    IKD::vector::Vec3Normalize(&N,&N);
    *pOut = V - N*((1+e)*IKD::vector::Vec3Dot(&N,&V));
}


// 壁との反射後の位置を算出
void GetRelectedPos( double Res_time, CIRCLE &circle, IKD::vector &RefV )
{
    // 衝突位置
    // 0.99は壁にめり込まないための補正
    circle.x = circle.Pre_x + circle.vx * (1-Res_time)*0.99f;
    circle.y = circle.Pre_y + circle.vy * (1-Res_time)*0.99f;
    // 反射ベクトル
    circle.vx = RefV.x;
    circle.vy = RefV.y;
    // 位置を補正
    circle.x += circle.vx * Res_time;
    circle.y += circle.vy * Res_time;
}


// 次の円の位置を算出
void GetNextCirclePos( CIRCLE &circle )
{
    IKD::vector RefV;    // 反射後の速度ベクトル
    double Res_time = 0.0f;    // 衝突後の移動可能時間

    // 重力を掛けて落とす
    circle.vy += g_Gravity/60;    // 1フレームで9.8/60(m/s)加速

    // 今の速度で位置を更新
    circle.Pre_x = circle.x;    // 前の位置を保存
    circle.Pre_y = circle.y;
    circle.x += circle.vx;        // 位置更新
    circle.y += circle.vy;

    // 壁との衝突をチェック
    // 左壁
    if(circle.x<0){
        // 反射後の速度ベクトルを取得
        GetRefrectVelo( &RefV, IKD::vector(1,0,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        // 残り時間算出
        Res_time = circle.x / circle.vx;
        // 反射後の位置を算出
        GetRelectedPos(    Res_time, circle, RefV );
    }
    // 右壁
    else if(circle.x>640){
        GetRefrectVelo( &RefV, IKD::vector(-1,0,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        Res_time = (circle.x-640) / circle.vx;
        GetRelectedPos(    Res_time, circle, RefV );
    }
    // 下壁
    else if(circle.y>480){
        GetRefrectVelo( &RefV, IKD::vector(0,-1,0), IKD::vector(circle.vx, circle.vy, 0), Wall_Ref);
        Res_time = (circle.y-480) / circle.vy;
        GetRelectedPos(    Res_time, circle, RefV );
    }
}


// 2円衝突処理
void CircleColProc( CIRCLE *c1, CIRCLE *c2 )
{
    double t=0;
    IKD::vector C1ColPos, C2ColPos, C1Velo, C2Velo;

    // 衝突している2円の衝突位置を検出
    if(!CalcParticleCollision(c1->r, c2->r,
        IKD::vector(c1->Pre_x, c1->Pre_y,0),
        IKD::vector(c1->x, c1->y, 0),
        IKD::vector(c2->Pre_x, c2->Pre_y, 0),
        IKD::vector(c2->x, c2->y, 0),
        t,
        C1ColPos,
        C2ColPos))
        return;    // 衝突していないようです

    // 衝突位置を前位置として保存
    c1->x = C1ColPos.x;
    c1->y = C1ColPos.y; 
    c2->x = C2ColPos.x;
    c2->y = C2ColPos.y; 
    c1->Pre_x = C1ColPos.x;
    c1->Pre_y = C1ColPos.y; 
    c2->Pre_x = C2ColPos.x;
    c2->Pre_y = C2ColPos.y; 

    // 衝突後の速度を算出
    if(!CalcParticleColliAfterPos(
        C1ColPos, IKD::vector(c1->vx, c1->vy, 0),
        C2ColPos, IKD::vector(c2->vx, c2->vy, 0),
        c1->w, c2->w,
        g_Circle_Ref, g_Circle_Ref,        // 球の反発係数
        t,
        C1ColPos, C1Velo,
        C2ColPos, C2Velo))
        return; // 何か失敗したようです

    // 衝突後位置に移動
    c1->vx = C1Velo.x;
    c1->vy = C1Velo.y;
    c2->vx = C2Velo.x;
    c2->vy = C2Velo.y;
    c1->x += c1->vx;
    c1->y += c1->vy;
    c2->x += c2->vx;
    c2->y += c2->vy;
}

// メイン関数
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPTSTR lpCmdLine, int nCmdShow){
    SetOutApplicationLogValidFlag(FALSE);
    ChangeWindowMode(TRUE);
    SetGraphMode(640, 480, 32);
    SetMainWindowText("2d collision");
    if(DxLib_Init() != 0){
        return -1;
    }

    int main_graphic_handle = MakeScreen(640, 480);

    srand(0x09485524);    // 乱数再初期化

    //WaitKey();

    int ff = sizeof(cell<CIRCLE>);

    /////////////////////////////
    // 円オブジェクトの初期化
    CIRCLE CAry[ g_CircleNum ];
    tree_object<CIRCLE> *spOFTAry[g_CircleNum];    // 円オブジェクトを包むOFTオブジェクト

    //////////////////////////////////////////
    // 円オブジェクトの初期位置・速度の設定
    //  遊ぶポイントです(^-^)
    //////////
    int cn;
    for(cn=0;cn<g_CircleNum;cn++)    // g_CircleNumだけ円を生成
    {
        CAry[cn].ID = cn;
        CAry[cn].r = 2 + 3*(double)rand()/RAND_MAX;    // 円の半径(2～5までランダム)
        CAry[cn].x = (double)cn/g_CircleNum*120+30*(double)rand()/RAND_MAX;        // 細長い初期位置
        CAry[cn].y = -400 + 700*(double)cn/g_CircleNum;    // 結構高いところから落とします
        CAry[cn].vx = 0.5;        // 初速（適当）
        CAry[cn].vy = 0;
        CAry[cn].scale = CAry[cn].r/32.0f;                // 画像が64×64なのでスケール値はこうなるんです！
        CAry[cn].w = CAry[cn].r*CAry[cn].r*CAry[cn].r;    // 質量は半径の3乗に比例とします
        // OFTに登録
        tree_object<CIRCLE> *p = new tree_object<CIRCLE>( cn );
        p->object = &CAry[cn];    // 登録
        spOFTAry[cn] = p;
    }


    ///////////////////////////////////
    // 線形4分木マネージャ
    //  空間範囲をX=-60～720; Y=-1200～520に設定
    //  円が飛び出さない範囲を指定すれば良いので
    //  アバウトです
    liner_for_tree_manager<CIRCLE> LTree;
    if(!LTree.init(g_PartitionLebel, -60.0f, -1200.0f, 720.0f, 520.0f))
    {
        return -1;
    }

    // ループ内一時変数
    std::uint32_t ColNum;                // 衝突判定回数
    collision_list<CIRCLE>* ColVect;    // 衝突対象リスト

    while(true){
        Sleep(16);
        if(ProcessMessage() != 0){
            break;
        }

        SetDrawScreen(main_graphic_handle);
        SetDrawScreen(DX_SCREEN_BACK);
        DrawBox(0, 0, 640, 480, GetColor(0xFF, 0xFF, 0xFF), TRUE);

        // 円の位置を更新してツリーに再登録
        static bool change = true;
        for( cn = 0; cn < g_CircleNum; cn++ )
        {
            int i = cn;
            // ↓テストのために逸脱再登録の順番を変えています
            if ( change == true ) {
                i = g_CircleNum - cn - 1;
            }
            CIRCLE *pTmp = spOFTAry[i]->object;
            GetNextCirclePos( *pTmp );    // 次の移動位置を仮決定
            spOFTAry[i]->remove();        // 一度リストから外れる
            // 再登録
            LTree.register_object( pTmp->x-pTmp->r, pTmp->y-pTmp->r, pTmp->x+pTmp->r, pTmp->y+pTmp->r, spOFTAry[i] );
        }
        change = !change;

        // 衝突対応リストを取得
        ColNum = LTree.get_all_collision_list( ColVect );

        // 衝突判定
        std::uint32_t c;
        ColNum/=2;    // 2で割るのはペアになっているので
        CIRCLE** pRoot = ColVect->root();
        for(c=0; c<ColNum; c++){
            double r2 = (pRoot[c*2]->r+pRoot[c*2+1]->r)*(pRoot[c*2]->r+pRoot[c*2+1]->r);
            double x = (pRoot[c*2]->x-pRoot[c*2+1]->x);
            double y = (pRoot[c*2]->y-pRoot[c*2+1]->y);
            if(r2 >= x*x + y*y )
            {
                // 2円衝突処理をする
                CircleColProc( pRoot[c*2], pRoot[c*2+1] );
            }
        }

        // 描画
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