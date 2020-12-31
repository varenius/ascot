#!/bin/sh
LD_LIBRARY_PATH=/data1/tobson/soft/ASCOT/qt5inst/qtbase/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
QT_PLUGIN_PATH=/data1/tobson/soft/ASCOT/qt5inst/qtbase/plugins${QT_PLUGIN_PATH:+:$QT_PLUGIN_PATH}
export QT_PLUGIN_PATH
exec /data1/tobson/soft/ASCOT/qt5inst/qtbase/bin/uic "$@"
