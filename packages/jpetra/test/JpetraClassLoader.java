/*
 * JpetraClassLoader.java
 *
 * Created on July 26, 2004, 4:03 PM
 */

package test;

import java.net.Socket;
import java.io.ObjectInputStream;
/**
 *
 * @author  jacross
 */
public class JpetraClassLoader extends ClassLoader {
    ObjectInputStream ois;
    
    public JpetraClassLoader() {
        // empty
    }
    
    public void setObjectInputStream(ObjectInputStream ois) {
        this.ois = ois;
    }
    
    public Class findClass(String name) {
        Class loadedClass = this.findLoadedClass(name);
        if (loadedClass != null) {
            return loadedClass;
        }
        
        byte[] b = loadClassData(name);
        return defineClass(name, b, 0, b.length);
    }
    
    private byte[] loadClassData(String name) {
        try {
            return (byte[]) ois.readObject();
        }
        catch (java.io.IOException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
        catch (java.lang.ClassNotFoundException e) {
            System.err.println(e.toString());
            System.exit(1);
        }
    
        return null;
    }
}
