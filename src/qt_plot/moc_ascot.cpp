/****************************************************************************
** Meta object code from reading C++ file 'ascot.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.12.8)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "ascot.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'ascot.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.12.8. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Ascot_t {
    QByteArrayData data[29];
    char stringdata0[427];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Ascot_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Ascot_t qt_meta_stringdata_Ascot = {
    {
QT_MOC_LITERAL(0, 0, 5), // "Ascot"
QT_MOC_LITERAL(1, 6, 8), // "_loadCfg"
QT_MOC_LITERAL(2, 15, 0), // ""
QT_MOC_LITERAL(3, 16, 16), // "_startLoadThread"
QT_MOC_LITERAL(4, 33, 3), // "row"
QT_MOC_LITERAL(5, 37, 6), // "column"
QT_MOC_LITERAL(6, 44, 19), // "_startAnalyseThread"
QT_MOC_LITERAL(7, 64, 21), // "_startResidualsThread"
QT_MOC_LITERAL(8, 86, 18), // "_startExportThread"
QT_MOC_LITERAL(9, 105, 19), // "_loadThreadFinished"
QT_MOC_LITERAL(10, 125, 22), // "_analyseThreadFinished"
QT_MOC_LITERAL(11, 148, 24), // "_residualsThreadFinished"
QT_MOC_LITERAL(12, 173, 21), // "_exportThreadFinished"
QT_MOC_LITERAL(13, 195, 20), // "_updateResidualTable"
QT_MOC_LITERAL(14, 216, 32), // "map<delaytype,vector<Residual> >"
QT_MOC_LITERAL(15, 249, 6), // "resids"
QT_MOC_LITERAL(16, 256, 19), // "_createResidualPlot"
QT_MOC_LITERAL(17, 276, 15), // "_appendLogInGUI"
QT_MOC_LITERAL(18, 292, 4), // "text"
QT_MOC_LITERAL(19, 297, 19), // "_getClickedTableRow"
QT_MOC_LITERAL(20, 317, 6), // "sender"
QT_MOC_LITERAL(21, 324, 17), // "_init_axis_layout"
QT_MOC_LITERAL(22, 342, 20), // "_setEliminationTable"
QT_MOC_LITERAL(23, 363, 13), // "ivg::Session&"
QT_MOC_LITERAL(24, 377, 7), // "session"
QT_MOC_LITERAL(25, 385, 8), // "_plotTrf"
QT_MOC_LITERAL(26, 394, 8), // "_plotCrf"
QT_MOC_LITERAL(27, 403, 17), // "closeTabWidgetTab"
QT_MOC_LITERAL(28, 421, 5) // "index"

    },
    "Ascot\0_loadCfg\0\0_startLoadThread\0row\0"
    "column\0_startAnalyseThread\0"
    "_startResidualsThread\0_startExportThread\0"
    "_loadThreadFinished\0_analyseThreadFinished\0"
    "_residualsThreadFinished\0_exportThreadFinished\0"
    "_updateResidualTable\0"
    "map<delaytype,vector<Residual> >\0"
    "resids\0_createResidualPlot\0_appendLogInGUI\0"
    "text\0_getClickedTableRow\0sender\0"
    "_init_axis_layout\0_setEliminationTable\0"
    "ivg::Session&\0session\0_plotTrf\0_plotCrf\0"
    "closeTabWidgetTab\0index"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Ascot[] = {

 // content:
       8,       // revision
       0,       // classname
       0,    0, // classinfo
      18,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  104,    2, 0x08 /* Private */,
       3,    2,  105,    2, 0x08 /* Private */,
       6,    0,  110,    2, 0x08 /* Private */,
       7,    0,  111,    2, 0x08 /* Private */,
       8,    0,  112,    2, 0x08 /* Private */,
       9,    0,  113,    2, 0x08 /* Private */,
      10,    0,  114,    2, 0x08 /* Private */,
      11,    0,  115,    2, 0x08 /* Private */,
      12,    0,  116,    2, 0x08 /* Private */,
      13,    1,  117,    2, 0x08 /* Private */,
      16,    0,  120,    2, 0x08 /* Private */,
      17,    1,  121,    2, 0x08 /* Private */,
      19,    1,  124,    2, 0x08 /* Private */,
      21,    0,  127,    2, 0x08 /* Private */,
      22,    1,  128,    2, 0x08 /* Private */,
      25,    1,  131,    2, 0x08 /* Private */,
      26,    1,  134,    2, 0x08 /* Private */,
      27,    1,  137,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    4,    5,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 14,   15,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString,   18,
    QMetaType::Int, QMetaType::QObjectStar,   20,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 23,   24,
    QMetaType::Void, 0x80000000 | 23,   24,
    QMetaType::Void, 0x80000000 | 23,   24,
    QMetaType::Void, QMetaType::Int,   28,

       0        // eod
};

void Ascot::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Ascot *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->_loadCfg(); break;
        case 1: _t->_startLoadThread((*reinterpret_cast< const int(*)>(_a[1])),(*reinterpret_cast< const int(*)>(_a[2]))); break;
        case 2: _t->_startAnalyseThread(); break;
        case 3: _t->_startResidualsThread(); break;
        case 4: _t->_startExportThread(); break;
        case 5: _t->_loadThreadFinished(); break;
        case 6: _t->_analyseThreadFinished(); break;
        case 7: _t->_residualsThreadFinished(); break;
        case 8: _t->_exportThreadFinished(); break;
        case 9: _t->_updateResidualTable((*reinterpret_cast< map<delaytype,vector<Residual> >(*)>(_a[1]))); break;
        case 10: _t->_createResidualPlot(); break;
        case 11: _t->_appendLogInGUI((*reinterpret_cast< const QString(*)>(_a[1]))); break;
        case 12: { int _r = _t->_getClickedTableRow((*reinterpret_cast< QObject*(*)>(_a[1])));
            if (_a[0]) *reinterpret_cast< int*>(_a[0]) = std::move(_r); }  break;
        case 13: _t->_init_axis_layout(); break;
        case 14: _t->_setEliminationTable((*reinterpret_cast< ivg::Session(*)>(_a[1]))); break;
        case 15: _t->_plotTrf((*reinterpret_cast< ivg::Session(*)>(_a[1]))); break;
        case 16: _t->_plotCrf((*reinterpret_cast< ivg::Session(*)>(_a[1]))); break;
        case 17: _t->closeTabWidgetTab((*reinterpret_cast< int(*)>(_a[1]))); break;
        default: ;
        }
    }
}

QT_INIT_METAOBJECT const QMetaObject Ascot::staticMetaObject = { {
    &QMainWindow::staticMetaObject,
    qt_meta_stringdata_Ascot.data,
    qt_meta_data_Ascot,
    qt_static_metacall,
    nullptr,
    nullptr
} };


const QMetaObject *Ascot::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Ascot::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Ascot.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int Ascot::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 18)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 18;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 18)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 18;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
