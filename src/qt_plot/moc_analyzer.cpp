/****************************************************************************
** Meta object code from reading C++ file 'analyzer.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.9.5)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "analyzer.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'analyzer.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.9.5. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_Analyzer_t {
    QByteArrayData data[31];
    char stringdata0[483];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_Analyzer_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_Analyzer_t qt_meta_stringdata_Analyzer = {
    {
QT_MOC_LITERAL(0, 0, 8), // "Analyzer"
QT_MOC_LITERAL(1, 9, 17), // "objectTypeChanged"
QT_MOC_LITERAL(2, 27, 0), // ""
QT_MOC_LITERAL(3, 28, 19), // "clearLoadedSessions"
QT_MOC_LITERAL(4, 48, 21), // "clearHighlightedItems"
QT_MOC_LITERAL(5, 70, 19), // "updateDiffsGroupBox"
QT_MOC_LITERAL(6, 90, 5), // "index"
QT_MOC_LITERAL(7, 96, 18), // "addObjectParameter"
QT_MOC_LITERAL(8, 115, 23), // "updateSelectedReference"
QT_MOC_LITERAL(9, 139, 3), // "row"
QT_MOC_LITERAL(10, 143, 6), // "column"
QT_MOC_LITERAL(11, 150, 19), // "selectedFileChanged"
QT_MOC_LITERAL(12, 170, 7), // "current"
QT_MOC_LITERAL(13, 178, 8), // "previous"
QT_MOC_LITERAL(14, 187, 18), // "plotReferenceFrame"
QT_MOC_LITERAL(15, 206, 13), // "loadEopSeries"
QT_MOC_LITERAL(16, 220, 22), // "performStationAnalysis"
QT_MOC_LITERAL(17, 243, 20), // "plotSelectedStations"
QT_MOC_LITERAL(18, 264, 16), // "plotSelectedEops"
QT_MOC_LITERAL(19, 281, 19), // "plotSelectedSources"
QT_MOC_LITERAL(20, 301, 23), // "loadLatestSelectedSinex"
QT_MOC_LITERAL(21, 325, 25), // "transformSelectedSessions"
QT_MOC_LITERAL(22, 351, 17), // "closeTabWidgetTab"
QT_MOC_LITERAL(23, 369, 18), // "analysisTabChanged"
QT_MOC_LITERAL(24, 388, 20), // "highlightCurrentPlot"
QT_MOC_LITERAL(25, 409, 13), // "stringCompare"
QT_MOC_LITERAL(26, 423, 6), // "string"
QT_MOC_LITERAL(27, 430, 1), // "l"
QT_MOC_LITERAL(28, 432, 1), // "r"
QT_MOC_LITERAL(29, 434, 25), // "_openDirectoryQFileDialog"
QT_MOC_LITERAL(30, 460, 22) // "_refreshFolderTreeview"

    },
    "Analyzer\0objectTypeChanged\0\0"
    "clearLoadedSessions\0clearHighlightedItems\0"
    "updateDiffsGroupBox\0index\0addObjectParameter\0"
    "updateSelectedReference\0row\0column\0"
    "selectedFileChanged\0current\0previous\0"
    "plotReferenceFrame\0loadEopSeries\0"
    "performStationAnalysis\0plotSelectedStations\0"
    "plotSelectedEops\0plotSelectedSources\0"
    "loadLatestSelectedSinex\0"
    "transformSelectedSessions\0closeTabWidgetTab\0"
    "analysisTabChanged\0highlightCurrentPlot\0"
    "stringCompare\0string\0l\0r\0"
    "_openDirectoryQFileDialog\0"
    "_refreshFolderTreeview"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Analyzer[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      21,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,  119,    2, 0x0a /* Public */,
       3,    0,  120,    2, 0x0a /* Public */,
       4,    0,  121,    2, 0x0a /* Public */,
       5,    1,  122,    2, 0x0a /* Public */,
       7,    0,  125,    2, 0x0a /* Public */,
       8,    2,  126,    2, 0x0a /* Public */,
      11,    2,  131,    2, 0x0a /* Public */,
      14,    0,  136,    2, 0x0a /* Public */,
      15,    0,  137,    2, 0x0a /* Public */,
      16,    0,  138,    2, 0x0a /* Public */,
      17,    0,  139,    2, 0x0a /* Public */,
      18,    0,  140,    2, 0x0a /* Public */,
      19,    0,  141,    2, 0x0a /* Public */,
      20,    0,  142,    2, 0x0a /* Public */,
      21,    0,  143,    2, 0x0a /* Public */,
      22,    1,  144,    2, 0x0a /* Public */,
      23,    1,  147,    2, 0x0a /* Public */,
      24,    0,  150,    2, 0x0a /* Public */,
      25,    2,  151,    2, 0x08 /* Private */,
      29,    0,  156,    2, 0x08 /* Private */,
      30,    0,  157,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::QModelIndex,    6,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int, QMetaType::Int,    9,   10,
    QMetaType::Void, QMetaType::QModelIndex, QMetaType::QModelIndex,   12,   13,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Int,    6,
    QMetaType::Void, QMetaType::Int,    6,
    QMetaType::Void,
    QMetaType::Bool, 0x80000000 | 26, 0x80000000 | 26,   27,   28,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void Analyzer::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Analyzer *_t = static_cast<Analyzer *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->objectTypeChanged(); break;
        case 1: _t->clearLoadedSessions(); break;
        case 2: _t->clearHighlightedItems(); break;
        case 3: _t->updateDiffsGroupBox((*reinterpret_cast< const QModelIndex(*)>(_a[1]))); break;
        case 4: _t->addObjectParameter(); break;
        case 5: _t->updateSelectedReference((*reinterpret_cast< const int(*)>(_a[1])),(*reinterpret_cast< const int(*)>(_a[2]))); break;
        case 6: _t->selectedFileChanged((*reinterpret_cast< const QModelIndex(*)>(_a[1])),(*reinterpret_cast< const QModelIndex(*)>(_a[2]))); break;
        case 7: _t->plotReferenceFrame(); break;
        case 8: _t->loadEopSeries(); break;
        case 9: _t->performStationAnalysis(); break;
        case 10: _t->plotSelectedStations(); break;
        case 11: _t->plotSelectedEops(); break;
        case 12: _t->plotSelectedSources(); break;
        case 13: _t->loadLatestSelectedSinex(); break;
        case 14: _t->transformSelectedSessions(); break;
        case 15: _t->closeTabWidgetTab((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 16: _t->analysisTabChanged((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 17: _t->highlightCurrentPlot(); break;
        case 18: { bool _r = _t->stringCompare((*reinterpret_cast< const string(*)>(_a[1])),(*reinterpret_cast< const string(*)>(_a[2])));
            if (_a[0]) *reinterpret_cast< bool*>(_a[0]) = std::move(_r); }  break;
        case 19: _t->_openDirectoryQFileDialog(); break;
        case 20: _t->_refreshFolderTreeview(); break;
        default: ;
        }
    }
}

const QMetaObject Analyzer::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_Analyzer.data,
      qt_meta_data_Analyzer,  qt_static_metacall, nullptr, nullptr}
};


const QMetaObject *Analyzer::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Analyzer::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Analyzer.stringdata0))
        return static_cast<void*>(this);
    return QMainWindow::qt_metacast(_clname);
}

int Analyzer::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 21)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 21;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 21)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 21;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
