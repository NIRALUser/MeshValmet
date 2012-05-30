/* $Id: TextWidget.cxx,v 1.1 2005/05/24 07:37:58 xushun Exp $ */


/*
 *
 *  Copyright (C) 2001-2004 EPFL (Swiss Federal Institute of Technology,
 *  Lausanne) This program is free software; you can redistribute it
 *  and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of
 *  the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 *  USA.
 *
 *  In addition, as a special exception, EPFL gives permission to link
 *  the code of this program with the Qt non-commercial edition library
 *  (or with modified versions of Qt non-commercial edition that use the
 *  same license as Qt non-commercial edition), and distribute linked
 *  combinations including the two.  You must obey the GNU General
 *  Public License in all respects for all of the code used other than
 *  Qt non-commercial edition.  If you modify this file, you may extend
 *  this exception to your version of the file, but you are not
 *  obligated to do so.  If you do not wish to do so, delete this
 *  exception statement from your version.
 *
 *  Authors : Nicolas Aspert, Diego Santa-Cruz and Davy Jacquet
 *
 *  Web site : http://mesh.epfl.ch
 *
 *  Reference :
 *   "MESH : Measuring Errors between Surfaces using the Hausdorff distance"
 *   in Proceedings of IEEE Intl. Conf. on Multimedia and Expo (ICME) 2002, 
 *   vol. I, pp. 705-708, available on http://mesh.epfl.ch
 *
 */








#include <TextWidget.h>
#include <qapplication.h>

TextWidget::TextWidget(QWidget *parent, const char *name)
  : QWidget(parent, name) {
  QFont font(QApplication::font());

  layout = new QGridLayout(this,2,1);
  view = new QTextView(this);

#if QT_VERSION >= 300 // Only in Qt >= 3.0
  view->setWordWrap(QTextEdit::NoWrap);
#endif

  view->setTextFormat(Qt::PlainText);
  font.setFamily("courier");
  font.setStyleHint(QFont::TypeWriter,QFont::PreferQuality);
  font.setFixedPitch(TRUE);
  font.setPointSize(9);
  view->setFont(font);
  layout->addWidget(view,0,0);
  butClose = new QPushButton("Close", this);
  layout->addWidget(butClose,1,0,Qt::AlignCenter);

  connect(butClose, SIGNAL(clicked()), this, SLOT(close()));

  setCaption("Mesh execution log");
}

TextWidget::~TextWidget()
{
  delete view;
  delete layout;
  delete butClose;
}

QSize TextWidget::sizeHint() const {
  return QSize(600,500);
}

void TextWidget::append(const QString &str) {
  // The append() method of QT's TextView is buggy. Use the recommended
  // workaround.
  view->setText(view->text()+str);
  qApp->processEvents(100);
}

void TextWidget_puts(void *out, const char *str) {
  TextWidget *tw;

  tw = (TextWidget*) out;
  tw->append(QString(str));
}
