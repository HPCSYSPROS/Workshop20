#ifndef PLUGINIOMGR_H
#define PLUGINIOMGR_H

#include "molfile_plugin.h"
#include "libmolfile_plugin.h"
//#include "hash.h"

class PluginIOMgr{
//public:
    //static const int MAX_PLUGINS;
private:
    //hast_t pluginhash;
    //int numPlugins;
    //very initial implementation, later this
    //should be changed
    molfile_plugin_t *psfpdbPlugin;

public:
    PluginIOMgr();
    ~PluginIOMgr();
    //molfile_plugin_t *getPlugin(const char *filetype);
    molfile_plugin_t *getPlugin() { return psfpdbPlugin; }
    void setPlugin(molfile_plugin_t *p) { psfpdbPlugin = p; }
};

#endif
