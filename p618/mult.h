// mult.h : main header file for the MULT application
//

#if !defined(AFX_MULT_H__B743E964_8F03_4E18_8FA1_124A68B3F922__INCLUDED_)
#define AFX_MULT_H__B743E964_8F03_4E18_8FA1_124A68B3F922__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// main symbols

/////////////////////////////////////////////////////////////////////////////
// CMultApp:
// See mult.cpp for the implementation of this class
//

class CMultApp : public CWinApp
{
public:
	CMultApp();

// Overrides
	// ClassWizard generated virtual function overrides
	//{{AFX_VIRTUAL(CMultApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

// Implementation

	//{{AFX_MSG(CMultApp)
		// NOTE - the ClassWizard will add and remove member functions here.
		//    DO NOT EDIT what you see in these blocks of generated code !
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ will insert additional declarations immediately before the previous line.

#endif // !defined(AFX_MULT_H__B743E964_8F03_4E18_8FA1_124A68B3F922__INCLUDED_)
