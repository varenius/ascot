/****************************************************************************
** Meta object code from reading C++ file 'lp_sked_gui.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "lp_sked_gui.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'lp_sked_gui.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_lps__GUI_t {
    QByteArrayData data[13];
    char stringdata0[143];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_lps__GUI_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_lps__GUI_t qt_meta_stringdata_lps__GUI = {
    {
QT_MOC_LITERAL(0, 0, 8), // "lps::GUI"
QT_MOC_LITERAL(1, 9, 13), // "foundSolution"
QT_MOC_LITERAL(2, 23, 0), // ""
QT_MOC_LITERAL(3, 24, 33), // "std::vector<lps::StationActiv..."
QT_MOC_LITERAL(4, 58, 8), // "activity"
QT_MOC_LITERAL(5, 67, 12), // "selectedCell"
QT_MOC_LITERAL(6, 80, 5), // "level"
QT_MOC_LITERAL(7, 86, 12), // "temporal_idx"
QT_MOC_LITERAL(8, 99, 9), // "lps::Path"
QT_MOC_LITERAL(9, 109, 4), // "rect"
QT_MOC_LITERAL(10, 114, 7), // "sta_idx"
QT_MOC_LITERAL(11, 122, 5), // "valid"
QT_MOC_LITERAL(12, 128, 14) // "print_skyplots"

    },
    "lps::GUI\0foundSolution\0\0"
    "std::vector<lps::StationActivity>\0"
    "activity\0selectedCell\0level\0temporal_idx\0"
    "lps::Path\0rect\0sta_idx\0valid\0"
    "print_skyplots"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_lps__GUI[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x0a /* Public */,
       5,    5,   32,    2, 0x0a /* Public */,
      12,    0,   43,    2, 0x0a /* Public */,

 // slots: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, 0x80000000 | 8, QMetaType::Int, QMetaType::Bool,    6,    7,    9,   10,   11,
    QMetaType::Void,

       0        // eod
};

void lps::GUI::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GUI *_t = static_cast<GUI *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->foundSolution((*reinterpret_cast< std::vector<lps::StationActivity>(*)>(_a[1]))); break;
        case 1: _t->selectedCell((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< lps::Path(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4])),(*reinterpret_cast< bool(*)>(_a[5]))); break;
        case 2: _t->print_skyplots(); break;
        default: ;
        }
    }
}

const QMetaObject lps::GUI::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_lps__GUI.data,
      qt_meta_data_lps__GUI,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *lps::GUI::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *lps::GUI::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_lps__GUI.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int lps::GUI::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
