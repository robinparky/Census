/* $Id: UsedClassFileFilter.java,v 1.11.2.2 2007/01/18 21:31:53 eric Exp $
 *
 * ProGuard -- shrinking, optimization, and obfuscation of Java class files.
 *
 * Copyright (c) 2002-2007 Eric Lafortune (eric@graphics.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
package proguard.shrink;

import proguard.classfile.*;
import proguard.classfile.visitor.*;


/**
 * This ClassFileVisitor delegates all its method calls to another
 * ClassFileVisitor, but only for ClassFile objects that are marked as used.
 *
 * @see UsageMarker
 *
 * @author Eric Lafortune
 */
public class UsedClassFileFilter
  implements ClassFileVisitor
{
    private UsageMarker      usageMarker;
    private ClassFileVisitor classFileVisitor;


    /**
     * Creates a new UsedClassFileFilter.
     * @param usageMarker      the usage marker that is used to mark the classes
     *                         and class members.
     * @param classFileVisitor the class file visitor to which the visiting
     *                         will be delegated.
     */
    public UsedClassFileFilter(UsageMarker      usageMarker,
                                ClassFileVisitor classFileVisitor)
    {
        this.usageMarker      = usageMarker;
        this.classFileVisitor = classFileVisitor;
    }


    // Implementations for ClassFileVisitor.

    public void visitProgramClassFile(ProgramClassFile programClassFile)
    {
        if (usageMarker.isUsed(programClassFile))
        {
            classFileVisitor.visitProgramClassFile(programClassFile);
        }
    }


    public void visitLibraryClassFile(LibraryClassFile libraryClassFile)
    {
        if (usageMarker.isUsed(libraryClassFile))
        {
            classFileVisitor.visitLibraryClassFile(libraryClassFile);
        }
    }
}
