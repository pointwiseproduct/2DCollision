// ColTrees.h (2007.5.13)

// �Փ˔���p�c���[�N���X�錾��

// �C�����
// 2007. 5. 13
// 2�̋��ʋ�Ԃ��Z�o���鎞�Ƀ��[�g���ԍ��������Z���Ă����o�O��
// �r���I�_���a�ɏC���B

// 2008. 12. 5
// GetMortonNumber���\�b�h�ŏ�����Ԃ����̐e��ԂɂȂ��Ă����o�O���C��


// 2009. 2. 7
// �X�}�[�g�|�C���^�ɂ��Ǘ�����߂��i�z�Q�Ɩ��̂��߁j
// Remove����Pre-Next�ڑ��̃o�O���C��


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
// ���ؓo�^�I�u�W�F�N�g(OFT)
//////////////////////////////////
template< class T>
class tree_object
{
public:
    cell<T> *m_pCell;            // �o�^���
    T *m_pObject;                // ����ΏۃI�u�W�F�N�g
    tree_object<T>* m_spPre;    // �O��tree_object�\����
    tree_object<T>* m_spNext;    // ����tree_object�\����
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
    // ���烊�X�g����O���
    bool Remove(){
        // ���łɈ�E���Ă��鎞�͏����I��
        if(!m_pCell)
            return false;

        // ������o�^���Ă����ԂɎ��g��ʒm
        if(!m_pCell->OnRemove( this ))
            return false;

        // ��E����
        // �O��̃I�u�W�F�N�g�����т���
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
// �Փ˃��X�g
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
    T**        root_;        // ���X�g���[�g�|�C���^
    size_t    pos_;
    size_t    mallocSize_;
};

/////////////////////////////////////
// ���`4���؋�ԊǗ��N���X
//////////////////////////////////
#define CLINER4TREEMANAGER_MAXLEVEL        9
template <class T>
class liner_for_tree_manager
{
protected:
    unsigned int m_uiDim;
    cell<T> **ppCellAry;    // ���`��ԃ|�C���^�z��
    unsigned int m_iPow[CLINER4TREEMANAGER_MAXLEVEL+1];    // �ׂ��搔�l�z��
    double m_fW;        // �̈��X����
    double m_fH;        // �̈��Y����
    double m_fLeft;    // �̈�̍����iX���ŏ��l�j
    double m_fTop;    // �̈�̏㑤�iY���ŏ��l�j
    double m_fUnit_W;        // �ŏ����x����Ԃ̕��P��
    double m_fUnit_H;        // �ŏ����x����Ԃ̍��P��
    std::uint32_t m_dwCellNum;        // ��Ԃ̐�
    unsigned int m_uiLevel;            // �ŉ��ʃ��x��
    collision_list< T > m_ColList;    // �Փ˃��X�g

public:
    // �R���X�g���N�^
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

        // �e���x���ł̋�Ԑ����Z�o
        int i;
        m_iPow[0] = 1;
        for(i=1;i<CLINER4TREEMANAGER_MAXLEVEL+1;i++)
            m_iPow[i] = m_iPow[i-1]*4;
    }

    // �f�X�g���N�^
    virtual ~liner_for_tree_manager()
    {
        std::uint32_t i;
        for(i=0; i<m_dwCellNum; i++){
            if(ppCellAry[i]!=nullptr)
                delete ppCellAry[i];
        }
        delete[] ppCellAry;
    }

    // ���`4���ؔz����\�z����
    bool Init( unsigned int Level, double left, double top, double right, double bottom )
    {
        // �ݒ�ō����x���ȏ�̋�Ԃ͍��Ȃ�
        if(Level>=CLINER4TREEMANAGER_MAXLEVEL)
            return false;

        // Level���x���i0��_�j�̔z��쐬
        m_dwCellNum = (m_iPow[Level+1]-1)/3;
        ppCellAry = new cell<T>*[m_dwCellNum];
        ZeroMemory( ppCellAry, sizeof(cell<T>*)*m_dwCellNum );

        // �̈��o�^
        m_fLeft = left;
        m_fTop = top;
        m_fW = right - left;
        m_fH = bottom - top;
        m_fUnit_W = m_fW/(1<<Level);
        m_fUnit_H = m_fH/(1<<Level);

        m_uiLevel = Level;

        return true;
    }

    // �I�u�W�F�N�g��o�^����
    bool Regist( double left, double top, double right, double bottom, tree_object<T> *spOFT )
    {
        // �I�u�W�F�N�g�̋��E�͈͂���o�^���[�g���ԍ����Z�o
        std::uint32_t Elem = GetMortonNumber( left, top, right, bottom );
        if(Elem < m_dwCellNum ){
            // ��Ԃ������ꍇ�͐V�K�쐬
            if( !ppCellAry[Elem] )
                CreateNewCell(Elem);
            return ppCellAry[Elem]->Push(spOFT);
        }
        return false;    // �o�^���s
    }

    // �Փ˔��胊�X�g���쐬����
    std::uint32_t GetAllCollisionList( collision_list<T>** colList )
    {
        // ���X�g�i�z��j�͕K�����������܂�
        m_ColList.resetPos();

        // ���[�g��Ԃ̑��݂��`�F�b�N
        if( ppCellAry[0] == nullptr )
            return 0;    // ��Ԃ����݂��Ă��Ȃ�

        // ���[�g��Ԃ�����
        list<T*> ColStac;
        GetCollisionList( 0, ColStac );

        *colList = &m_ColList;

        return (std::uint32_t)m_ColList.getSize();
    }



protected:

    // ��ԓ��ŏՓ˃��X�g���쐬����
    bool GetCollisionList( std::uint32_t Elem, list<T*> &ColStac )
    {
        list<T*>::iterator it;
        // �@ ��ԓ��̃I�u�W�F�N�g���m�̏Փ˃��X�g�쐬
        tree_object<T>* spOFT1 = ppCellAry[Elem]->GetFirstObj();
        while( spOFT1 != 0 )
        {
            tree_object<T>* spOFT2 = spOFT1->m_spNext;
            while( spOFT2 != 0 ){
                m_ColList.write( spOFT1->m_pObject, spOFT2->m_pObject );
                spOFT2 = spOFT2->m_spNext;
            }
            // �A �Փ˃X�^�b�N�Ƃ̏Փ˃��X�g�쐬
            for( it = ColStac.begin(); it != ColStac.end(); it++ ){
                m_ColList.write( spOFT1->m_pObject, *it );
            }
            spOFT1 = spOFT1->m_spNext;
        }

        bool ChildFlag = false;
        // �B �q��ԂɈړ�
        std::uint32_t ObjNum = 0;
        std::uint32_t i, NextElem;
        for(i=0; i<4; i++){ 
            NextElem = Elem*4+1+i;
            if( NextElem<m_dwCellNum && ppCellAry[Elem*4+1+i] ){
                if(!ChildFlag){
                    // �C �o�^�I�u�W�F�N�g���X�^�b�N�ɒǉ�
                    spOFT1 = ppCellAry[Elem]->GetFirstObj();
                    while( spOFT1 != 0 ){
                        ColStac.push_back( spOFT1->m_pObject );
                        ObjNum++;
                        spOFT1 = spOFT1->m_spNext;
                    }
                }
                ChildFlag = true;
                GetCollisionList( Elem*4+1+i, ColStac );    // �q��Ԃ�
            }
        }

        // �D �X�^�b�N����I�u�W�F�N�g���O��
        if( ChildFlag ){
            for(i=0; i<ObjNum; i++)
                ColStac.pop_back();
        }

        return true;
    }


    // ��Ԃ𐶐�
    bool CreateNewCell( std::uint32_t Elem )
    {
        // �����̗v�f�ԍ�
        while( !ppCellAry[Elem] )
        {
            // �w��̗v�f�ԍ��ɋ�Ԃ�V�K�쐬
            ppCellAry[Elem] = new cell<T>;

            // �e��ԂɃW�����v
            Elem = (Elem-1)>>2;
            if(Elem>=m_dwCellNum) break;
        }
        return true;
    }

    // ���W�����Ԕԍ����Z�o
    std::uint32_t GetMortonNumber( double left, double top, double right, double bottom )
    {
        // �ŏ����x���ɂ�����e���ʒu���Z�o
        std::uint32_t LT = GetPointElem(left, top);
        std::uint32_t RB = GetPointElem(right, bottom );

        // ��Ԕԍ��̔r���I�_���a����
        // �������x�����Z�o
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

    // �r�b�g�����֐�
    std::uint32_t BitSeparate32( std::uint32_t n )
    {
        n = (n|(n<<8)) & 0x00ff00ff;
        n = (n|(n<<4)) & 0x0f0f0f0f;
        n = (n|(n<<2)) & 0x33333333;
        return (n|(n<<1)) & 0x55555555;
    }

    // 2D���[�g����Ԕԍ��Z�o�֐�
    std::uint16_t Get2DMortonNumber( std::uint16_t x, std::uint16_t y )
    {
        return (std::uint16_t)(BitSeparate32(x) | (BitSeparate32(y)<<1));
    }

    // ���W�����`4���ؗv�f�ԍ��ϊ��֐�
    std::uint32_t GetPointElem( double pos_x, double pos_y )
    {
        return Get2DMortonNumber( (std::uint16_t)((pos_x-m_fLeft)/m_fUnit_W), (std::uint16_t)((pos_y-m_fTop)/m_fUnit_H) );
    }
};


/////////////////////////////////////
// ��ԃN���X
//////////////////////////////////
template <class T>
class cell
{
protected:
    tree_object<T>* m_spLatest;    // �ŐVOFL�ւ̃X�}�[�g�|�C���^

public:
    // �R���X�g���N�^
    cell() : m_spLatest( 0 )
    {
    }

    // �f�X�g���N�^
    virtual ~cell()
    {
        if( m_spLatest != nullptr )
            ResetLink( m_spLatest );
    }

    // �����N��S�ă��Z�b�g����
    void ResetLink( tree_object<T> *spOFT )
    {
        if( spOFT->m_spNext != 0 ) {
            ResetLink( spOFT->m_spNext );
        }
        spOFT = 0;
    }

    // OFT���v�b�V��
    bool Push( tree_object<T> *spOFT )
    {
        if( spOFT == 0 ) return false;    // �����I�u�W�F�N�g�͓o�^���Ȃ�
        if( spOFT->m_pCell == this ) return false;    // 2�d�o�^�`�F�b�N
        if( m_spLatest == 0 ){
            m_spLatest = spOFT;    // ��ԂɐV�K�o�^
        }
        else
        {
            // �ŐVOFT�I�u�W�F�N�g���X�V
            spOFT->m_spNext = m_spLatest;
            m_spLatest->m_spPre = spOFT;
            m_spLatest = spOFT;
        }
        spOFT->m_pCell = this;    // ��Ԃ�o�^
        return true;
    }

    tree_object<T> *GetFirstObj()
    {
        return m_spLatest;
    }

    // �폜�����I�u�W�F�N�g���`�F�b�N
    bool OnRemove( tree_object<T> *pRemoveObj )
    {
        if( m_spLatest == pRemoveObj ){
            // ���̃I�u�W�F�N�g�ɑ}���ւ�
            if( m_spLatest != nullptr )
                m_spLatest = m_spLatest->m_spNext;
        }
        return true;
    }
};


}  // end namespace IKD