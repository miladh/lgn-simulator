TEMPLATE = subdirs
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += ordered

SUBDIRS += src apps tests \
    apps/firingSynchrony

OTHER_FILES += .qmake.conf .gitignore README.md


