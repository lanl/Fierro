#pragma once
#ifndef YAML_POST_PROCESSING_H
#define YAML_POST_PROCESSING_H

namespace Yaml {
    /**
     * Base class tags for enabling certain functionality.
     * 
     * DerivedFields:
     *  Will invoke obj.derive() immediately after serialization.
     * 
     * ValidatedYaml:
     *  Will invoke obj.validate() after serialization and after deriving any fields.
     * 
    */
    struct DerivedFields {
        virtual void derive() { }
    };
    struct ValidatedYaml {
        virtual void validate() { }
    };
    
    /**
     * Invoke derive.
     * Ensures only the desired class's derive 
     * is invoked, and not derived classes' derive.
    */
    template<typename T>
    void derive(T& v) {
        v.T::derive();
    }

    /**
     * Invoke validate.
     * Ensures only the desired class's validate 
     * is invoked, and not derived classes' validate.
    */
    template<typename T>
    void validate(T& v) {
        v.T::validate();
    }
}
#endif