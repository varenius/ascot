/****************************************************************************
** Meta object code from reading C++ file 'lp_sked_worker.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "lp_sked_worker.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'lp_sked_worker.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_lps__Worker_t {
    QByteArrayData data[13];
    char stringdata0[139];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_lps__Worker_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_lps__Worker_t qt_meta_stringdata_lps__Worker = {
    {
QT_MOC_LITERAL(0, 0, 11), // "lps::Worker"
QT_MOC_LITERAL(1, 12, 13), // "foundSolution"
QT_MOC_LITERAL(2, 26, 0), // ""
QT_MOC_LITERAL(3, 27, 33), // "std::vector<lps::StationActiv..."
QT_MOC_LITERAL(4, 61, 8), // "activity"
QT_MOC_LITERAL(5, 70, 12), // "selectedCell"
QT_MOC_LITERAL(6, 83, 5), // "level"
QT_MOC_LITERAL(7, 89, 12), // "temporal_idx"
QT_MOC_LITERAL(8, 102, 9), // "lps::Path"
QT_MOC_LITERAL(9, 112, 4), // "rect"
QT_MOC_LITERAL(10, 117, 7), // "station"
QT_MOC_LITERAL(11, 125, 5), // "valid"
QT_MOC_LITERAL(12, 131, 7) // "process"

    },
    "lps::Worker\0foundSolution\0\0"
    "std::vector<lps::StationActivity>\0"
    "activity\0selectedCell\0level\0temporal_idx\0"
    "lps::Path\0rect\0station\0valid\0process"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_lps__Worker[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       2,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   29,    2, 0x06 /* Public */,
       5,    5,   32,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
      12,    0,   43,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    4,
    QMetaType::Void, QMetaType::Int, QMetaType::Int, 0x80000000 | 8, QMetaType::Int, QMetaType::Bool,    6,    7,    9,   10,   11,

 // slots: parameters
    QMetaType::Void,

       0        // eod
};

void lps::Worker::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Worker *_t = static_cast<Worker *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->foundSolution((*reinterpret_cast< std::vector<lps::StationActivity>(*)>(_a[1]))); break;
        case 1: _t->selectedCell((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2])),(*reinterpret_cast< lps::Path(*)>(_a[3])),(*reinterpret_cast< int(*)>(_a[4])),(*reinterpret_cast< bool(*)>(_a[5]))); break;
        case 2: _t->process(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            typedef void (Worker::*_t)(std::vector<lps::StationActivity> );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Worker::foundSolution)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (Worker::*_t)(int , int , lps::Path , int , bool );
            if (*reinterpret_cast<_t *>(_a[1]) == static_cast<_t>(&Worker::selectedCell)) {
                *result = 1;
                return;
            }
        }
    }
}

const QMetaObject lps::Worker::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_lps__Worker.data,
      qt_meta_data_lps__Worker,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *lps::Worker::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *lps::Worker::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_lps__Worker.stringdata0))
        return static_cast<void*>(this);
    return QObject::qt_metacast(_clname);
}

int lps::Worker::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
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

// SIGNAL 0
void lps::Worker::foundSolution(std::vector<lps::StationActivity> _t1)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}

// SIGNAL 1
void lps::Worker::selectedCell(int _t1, int _t2, lps::Path _t3, int _t4, bool _t5)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)), const_cast<void*>(reinterpret_cast<const void*>(&_t4)), const_cast<void*>(reinterpret_cast<const void*>(&_t5)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
