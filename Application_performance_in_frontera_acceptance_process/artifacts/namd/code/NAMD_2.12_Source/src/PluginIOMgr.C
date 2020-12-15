#include <stdio.h>
#include "PluginIOMgr.h"

//static const int PluginIOMgr::MAX_PLUGINS = 200;

int register_cb(void *v, vmdplugin_t *p){
/*    const char *key = p->name;
    if (num_plugins >= MAX_PLUGINS) {
        fprintf(stderr, "Exceeded maximum allowed number of plugins; recompile. :(\n");
        return NAMDPLUGIN_ERROR;
    }
    if (hash_insert(&pluginhash, key, num_plugins) != HASH_FAIL) {
        fprintf(stderr, "Multiple plugins for file type '%s' found!", key);
        return NAMDPLUGIN_ERROR;
    }
    plugins[num_plugins++] = (molfile_plugin_t *)p;
    return NAMDPLUGIN_SUCCESS;
*/
    PluginIOMgr *me = (PluginIOMgr *)v;
    me->setPlugin((molfile_plugin_t *)p);
    return 0;
}

PluginIOMgr::PluginIOMgr(){
    molfile_jsplugin_init();
    molfile_jsplugin_register(this, register_cb);
}

PluginIOMgr::~PluginIOMgr(){
    molfile_jsplugin_fini();
}

/*molfile_plugin_t *PulginIOMgr::getPlugin(const char *filetype){
    int id;
    if ((id = hash_lookup(&pluginhash, filetype)) == HASH_FAIL) {
        fprintf(stderr, "No plugin found for filetype '%s'\n", filetype);
        return NULL;
    }

    return plugins[id];
}*/
