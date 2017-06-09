// CnvSSwin.h : main header file for the CNVSSWIN DLL
//

#if !defined(AFX_CNVSSWIN_H__552A1313_6694_11D1_9F23_0000F823B458__INCLUDED_)
#define AFX_CNVSSWIN_H__552A1313_6694_11D1_9F23_0000F823B458__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#ifndef __AFXWIN_H__
#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols

/////////////////////////////////////////////////////////////////////////////
// CCnvSSwinApp
// See CnvSSwin.cpp for the implementation of this class
//

class CCnvSSwinApp : public CWinApp
{
public:
    CCnvSSwinApp();

// Overrides
    // ClassWizard generated virtual function overrides
    //{{AFX_VIRTUAL(CCnvSSwinApp)
    //}}AFX_VIRTUAL

    //{{AFX_MSG(CCnvSSwinApp)
    // NOTE - the ClassWizard will add and remove member functions here.
    //    DO NOT EDIT what you see in these blocks of generated code !
    //}}AFX_MSG
    DECLARE_MESSAGE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Developer Studio will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_CNVSSWIN_H__552A1313_6694_11D1_9F23_0000F823B458__INCLUDED_)
