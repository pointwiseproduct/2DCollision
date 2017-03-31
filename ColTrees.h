// ColTrees.h (2007.5.13)

// 衝突判定用ツリークラス宣言部

// 修正情報
// 2007. 5. 13
// 2つの共通空間を算出する時にモートン番号を引き算していたバグを
// 排他的論理和に修正。

// 2008. 12. 5
// GetMortonNumberメソッドで所属空間が一つ上の親空間になっていたバグを修正


// 2009. 2. 7
// スマートポインタによる管理をやめた（循環参照問題のため）
// Remove時のPre-Next接続のバグを修正


#pragma once

#include <set>
#include <vector>
#include <list>
#include <map>
#include <cstdint>

using namespace std;


namespace IKD
{


template <class T>
class cell;
/////////////////////////////////////
// 分木登録オブジェクト(OFT)
//////////////////////////////////
template< class T>
class tree_object
{
public:
    cell<T> *m_pCell;            // 登録空間
    T *m_pObject;                // 判定対象オブジェクト
    tree_object<T>* m_spPre;    // 前のtree_object構造体
    tree_object<T>* m_spNext;    // 次のtree_object構造体
    int m_id;

public:
    tree_object( int id = 0 ){
        m_id = id;
        m_pCell = nullptr;
        m_pObject = nullptr;
        m_spPre = nullptr;
        m_spNext = nullptr;
    }

    virtual ~tree_object()
    {
    }

public:
    // 自らリストから外れる
    bool Remove(){
        // すでに逸脱している時は処理終了
        if(!m_pCell)
            return false;

        // 自分を登録している空間に自身を通知
        if(!m_pCell->OnRemove( this ))
            return false;

        // 逸脱処理
        // 前後のオブジェクトを結びつける
        if(m_spPre != 0)
        {
            m_spPre->m_spNext = m_spNext;
        }
        if(m_spNext != 0)
        {
            m_spNext->m_spPre = m_spPre;
        }
        m_spPre = 0;
        m_spNext = 0;
        m_pCell = nullptr;
        return true;
    }

};


/////////////////////////////////////
// 衝突リスト
//////////////////////////////////
#define COLLISIONLIST_REALLOCSIZE    100
template < class T >
class collision_list {
public:
    collision_list() : root_( 0 ), pos_( 0 ), mallocSize_( 0 ) {
        root_ = (T**)malloc( 0 );
    }
    ~collision_list() {
        reflesh();
    }

    size_t getSize() {
        return pos_;
    }

    T** getRootPtr() {
        return root_;
    }

    void resetPos() {
        pos_ = 0;
    }

    void write( T* obj1, T*obj2 ) {
        if ( pos_ >= mallocSize_ ) {
            root_ = (T**)realloc( root_, sizeof (T*) * ( mallocSize_ + COLLISIONLIST_REALLOCSIZE ) );
            mallocSize_ += COLLISIONLIST_REALLOCSIZE;
        }
        root_[ pos_++ ] = obj1;
        root_[ pos_++ ] = obj2;
    }

    void reflesh() {
        if ( root_ ) {
            free( root_ );
        }
    }

private:
    T**        root_;        // リストルートポインタ
    size_t    pos_;
    size_t    mallocSize_;
};

/////////////////////////////////////
// 線形4分木空間管理クラス
//////////////////////////////////
#define CLINER4TREEMANAGER_MAXLEVEL        9
template <class T>
class liner_for_tree_manager
{
protected:
    unsigned int m_uiDim;
    cell<T> **ppCellAry;    // 線形空間ポインタ配列
    unsigned int m_iPow[CLINER4TREEMANAGER_MAXLEVEL+1];    // べき乗数値配列
    double m_fW;        // 領域のX軸幅
    double m_fH;        // 領域のY軸幅
    double m_fLeft;    // 領域の左側（X軸最小値）
    double m_fTop;    // 領域の上側（Y軸最小値）
    double m_fUnit_W;        // 最小レベル空間の幅単位
    double m_fUnit_H;        // 最小レベル空間の高単位
    std::uint32_t m_dwCellNum;        // 空間の数
    unsigned int m_uiLevel;            // 最下位レベル
    collision_list< T > m_ColList;    // 衝突リスト

public:
    // コンストラクタ
    liner_for_tree_manager()
    {
        m_uiLevel = 0;
        m_fW = 0.0;
        m_fH = 0.0;
        m_fLeft = 0.0;
        m_fTop = 0.0;
        m_fUnit_W = 0.0;
        m_fUnit_H = 0.0;
        m_dwCellNum = 0;
        ppCellAry = nullptr;
        m_uiDim = 0;

        // 各レベルでの空間数を算出
        int i;
        m_iPow[0] = 1;
        for(i=1;i<CLINER4TREEMANAGER_MAXLEVEL+1;i++)
            m_iPow[i] = m_iPow[i-1]*4;
    }

    // デストラクタ
    virtual ~liner_for_tree_manager()
    {
        std::uint32_t i;
        for(i=0; i<m_dwCellNum; i++){
            if(ppCellAry[i]!=nullptr)
                delete ppCellAry[i];
        }
        delete[] ppCellAry;
    }

    // 線形4分木配列を構築する
    bool Init( unsigned int Level, double left, double top, double right, double bottom )
    {
        // 設定最高レベル以上の空間は作れない
        if(Level>=CLINER4TREEMANAGER_MAXLEVEL)
            return false;

        // Levelレベル（0基点）の配列作成
        m_dwCellNum = (m_iPow[Level+1]-1)/3;
        ppCellAry = new cell<T>*[m_dwCellNum];
        ZeroMemory( ppCellAry, sizeof(cell<T>*)*m_dwCellNum );

        // 領域を登録
        m_fLeft = left;
        m_fTop = top;
        m_fW = right - left;
        m_fH = bottom - top;
        m_fUnit_W = m_fW/(1<<Level);
        m_fUnit_H = m_fH/(1<<Level);

        m_uiLevel = Level;

        return true;
    }

    // オブジェクトを登録する
    bool Regist( double left, double top, double right, double bottom, tree_object<T> *spOFT )
    {
        // オブジェクトの境界範囲から登録モートン番号を算出
        std::uint32_t Elem = GetMortonNumber( left, top, right, bottom );
        if(Elem < m_dwCellNum ){
            // 空間が無い場合は新規作成
            if( !ppCellAry[Elem] )
                CreateNewCell(Elem);
            return ppCellAry[Elem]->Push(spOFT);
        }
        return false;    // 登録失敗
    }

    // 衝突判定リストを作成する
    std::uint32_t GetAllCollisionList( collision_list<T>** colList )
    {
        // リスト（配列）は必ず初期化します
        m_ColList.resetPos();

        // ルート空間の存在をチェック
        if( ppCellAry[0] == nullptr )
            return 0;    // 空間が存在していない

        // ルート空間を処理
        list<T*> ColStac;
        GetCollisionList( 0, ColStac );

        *colList = &m_ColList;

        return (std::uint32_t)m_ColList.getSize();
    }



protected:

    // 空間内で衝突リストを作成する
    bool GetCollisionList( std::uint32_t Elem, list<T*> &ColStac )
    {
        list<T*>::iterator it;
        // ① 空間内のオブジェクト同士の衝突リスト作成
        tree_object<T>* spOFT1 = ppCellAry[Elem]->GetFirstObj();
        while( spOFT1 != 0 )
        {
            tree_object<T>* spOFT2 = spOFT1->m_spNext;
            while( spOFT2 != 0 ){
                m_ColList.write( spOFT1->m_pObject, spOFT2->m_pObject );
                spOFT2 = spOFT2->m_spNext;
            }
            // ② 衝突スタックとの衝突リスト作成
            for( it = ColStac.begin(); it != ColStac.end(); it++ ){
                m_ColList.write( spOFT1->m_pObject, *it );
            }
            spOFT1 = spOFT1->m_spNext;
        }

        bool ChildFlag = false;
        // ③ 子空間に移動
        std::uint32_t ObjNum = 0;
        std::uint32_t i, NextElem;
        for(i=0; i<4; i++){ 
            NextElem = Elem*4+1+i;
            if( NextElem<m_dwCellNum && ppCellAry[Elem*4+1+i] ){
                if(!ChildFlag){
                    // ④ 登録オブジェクトをスタックに追加
                    spOFT1 = ppCellAry[Elem]->GetFirstObj();
                    while( spOFT1 != 0 ){
                        ColStac.push_back( spOFT1->m_pObject );
                        ObjNum++;
                        spOFT1 = spOFT1->m_spNext;
                    }
                }
                ChildFlag = true;
                GetCollisionList( Elem*4+1+i, ColStac );    // 子空間へ
            }
        }

        // ⑤ スタックからオブジェクトを外す
        if( ChildFlag ){
            for(i=0; i<ObjNum; i++)
                ColStac.pop_back();
        }

        return true;
    }


    // 空間を生成
    bool CreateNewCell( std::uint32_t Elem )
    {
        // 引数の要素番号
        while( !ppCellAry[Elem] )
        {
            // 指定の要素番号に空間を新規作成
            ppCellAry[Elem] = new cell<T>;

            // 親空間にジャンプ
            Elem = (Elem-1)>>2;
            if(Elem>=m_dwCellNum) break;
        }
        return true;
    }

    // 座標から空間番号を算出
    std::uint32_t GetMortonNumber( double left, double top, double right, double bottom )
    {
        // 最小レベルにおける各軸位置を算出
        std::uint32_t LT = GetPointElem(left, top);
        std::uint32_t RB = GetPointElem(right, bottom );

        // 空間番号の排他的論理和から
        // 所属レベルを算出
        std::uint32_t Def = RB ^ LT;
        unsigned int HiLevel = 0;
        unsigned int i;
        for(i=0; i<m_uiLevel; i++)
        {
            std::uint32_t Check = (Def>>(i*2)) & 0x3;
            if( Check != 0 )
                HiLevel = i+1;
        }
        std::uint32_t SpaceNum = RB>>(HiLevel*2);
        std::uint32_t AddNum = (m_iPow[m_uiLevel-HiLevel]-1)/3;
        SpaceNum += AddNum;

        if(SpaceNum > m_dwCellNum)
            return 0xffffffff;

        return SpaceNum;
    }

    // ビット分割関数
    std::uint32_t BitSeparate32( std::uint32_t n )
    {
        n = (n|(n<<8)) & 0x00ff00ff;
        n = (n|(n<<4)) & 0x0f0f0f0f;
        n = (n|(n<<2)) & 0x33333333;
        return (n|(n<<1)) & 0x55555555;
    }

    // 2Dモートン空間番号算出関数
    std::uint16_t Get2DMortonNumber( std::uint16_t x, std::uint16_t y )
    {
        return (std::uint16_t)(BitSeparate32(x) | (BitSeparate32(y)<<1));
    }

    // 座標→線形4分木要素番号変換関数
    std::uint32_t GetPointElem( double pos_x, double pos_y )
    {
        return Get2DMortonNumber( (std::uint16_t)((pos_x-m_fLeft)/m_fUnit_W), (std::uint16_t)((pos_y-m_fTop)/m_fUnit_H) );
    }
};


/////////////////////////////////////
// 空間クラス
//////////////////////////////////
template <class T>
class cell
{
protected:
    tree_object<T>* m_spLatest;    // 最新OFLへのスマートポインタ

public:
    // コンストラクタ
    cell() : m_spLatest( 0 )
    {
    }

    // デストラクタ
    virtual ~cell()
    {
        if( m_spLatest != nullptr )
            ResetLink( m_spLatest );
    }

    // リンクを全てリセットする
    void ResetLink( tree_object<T> *spOFT )
    {
        if( spOFT->m_spNext != 0 ) {
            ResetLink( spOFT->m_spNext );
        }
        spOFT = 0;
    }

    // OFTをプッシュ
    bool Push( tree_object<T> *spOFT )
    {
        if( spOFT == 0 ) return false;    // 無効オブジェクトは登録しない
        if( spOFT->m_pCell == this ) return false;    // 2重登録チェック
        if( m_spLatest == 0 ){
            m_spLatest = spOFT;    // 空間に新規登録
        }
        else
        {
            // 最新OFTオブジェクトを更新
            spOFT->m_spNext = m_spLatest;
            m_spLatest->m_spPre = spOFT;
            m_spLatest = spOFT;
        }
        spOFT->m_pCell = this;    // 空間を登録
        return true;
    }

    tree_object<T> *GetFirstObj()
    {
        return m_spLatest;
    }

    // 削除されるオブジェクトをチェック
    bool OnRemove( tree_object<T> *pRemoveObj )
    {
        if( m_spLatest == pRemoveObj ){
            // 次のオブジェクトに挿げ替え
            if( m_spLatest != nullptr )
                m_spLatest = m_spLatest->m_spNext;
        }
        return true;
    }
};


}  // end namespace IKD