---
title: 监测程序运行状态并发送邮件
tags: Code Shell Linux
layout: article
license: true
toc: true
key: a20210519
pageview: true
header:
  theme: dark
  background: 'linear-gradient(135deg, rgb(34, 139, 87), rgb(139, 34, 139))'
article_header:
  type: overlay
  theme: dark
  background_color: '#123'
  background_image: false
aside:
    toc: true
sitemap: true
mathjax: true
author: YuXuan
show_author_profile: true
---
这篇博客整理一个脚本,用来在服务器上给自己发邮件,主要结合程序执行监测,在程序执行完毕之后给自己发邮件提醒,可以让自己及时查看结果并进行下一步的计算过着执行别的任务,做这个事情最主要的原因是想节省时间,同时也有自己懒的原因
{:.info}
<!--more-->

# 发送文字
在前面[监测程序运行的Shell脚本](https://yxli8023.github.io/2021/05/15/Shell-Monitor.html)这篇博客中,虽然提供了监测程序运行的脚本,但是任务是否执行完毕还是需要自己登录服务器进行查看,这里就想整理一个脚本,在程序执行完毕之后通过邮件方式进行提醒,这样可以让自己及时的查看计算结果并开始下一步的计划,极大的方便了自己,同时也可以节省不少的时间,专心去做别的事情.
```python
import smtplib
from email.mime.text import MIMEText
import sys,os
 
#QQ邮箱提供的SMTP服务器
mail_host = 'smtp.qq.com'
#服务器端口
port = 465
send_by = '********@qq.com' # QQ邮箱账号
password = 'asdfghjkl' # 开启QQ邮箱STMP后获得的一串符号
send_to = 'abcd@mail.com' # 目标邮箱
def send_email(title,content,):
    #创建了MIMEText类，相当于在写邮件内容，是plain类型
      message = MIMEText(content,'plain','utf-8')
      message["From"] = send_by
      message['To'] = send_to
      message['Subject'] = title
      try:
          #注意第三个参数，设置了转码的格式(我不设的时候会报解码错误)
          smpt = smtplib.SMTP_SSL(mail_host, port, 'utf-8')
          smpt.login(send_by,password)
          smpt.sendmail(send_by, send_to,message.as_string())
          print("发送成功")
      except:
          print("发送失败")
def main():
    title = sys.argv[1]  # 接受第一个输入参数
    content = sys.argv[2]  # 接受第二个输入参数
    send_email(title,content)
if __name__ == "__main__":
    main()
```
这里只需要配置好自己要发送邮件的qq邮箱与接受邮件消息的邮箱即可,qq邮箱需要去开启STMP服务,这个可以自行百度解决.配置号邮箱之后,就只需要确定输入参数邮件标题`title`与邮件内容`content`.这里使用了`sys.argv`这个函数,在使用这个python脚本的时候需要输入参数
```shell
python3 filename.py "title" "content"
```
通过这样的执行之后,邮件标题即为**title**,邮件内容为**content**.

# 实例
```shell
#!/bin/sh  
#============ get the file name ===========  
Folder_A=$(pwd) 
for file_a in ${Folder_A}/*.f90
	do 
	if [ -n "$file_a" ];then
		temp_file=`basename $file_a  .f90` 
		ifort -mkl -O3 -CB $file_a -o $temp_file.out 
		if [ -e $temp_file.out ];then
			./$temp_file.out &
		fi
	fi
done
python3 pathToMail/mail.py "Code Complier" "Your Code have been finished"
```
通过这个shell脚本,当文件夹中所有的fortran程序都编译执行完毕之后,就会发送一封邮件到我的邮箱提醒我任务进度.同样也可将这个发送邮件的python脚本和我前面监测程序运行的shell脚本结合,提醒程序执行进度,nice!!!!!

# 发送文字和图片
```python
# coding:utf-8
# Ref:https://blog.csdn.net/daimashiren/article/details/98470427
import smtplib
import os,sys
from email.mime.text import MIMEText # 文字
from email.mime.multipart import MIMEMultipart 
from email.header import Header
from email import encoders
from email.mime.base import MIMEBase
from email.utils import parseaddr, formataddr
import base64
from email.mime.image import MIMEImage # 图片
import traceback
 
 
class SendMail(object):
    def __init__(self, title=None, content=None,
                 file=None, image=None, ssl_port=465):
 
        '''
               :param username: 用户名
               :param passwd: 密码(QQ邮箱SMTP开启后的一串字符)
               :param recv: 收件人，可以同时给多个人发送,此时需要使用list ['user1@qq.com','user2@qq.com]
               :param title: 邮件标题
               :param content: 邮件正文
               :param image: 图片路径，绝对路径，默认为无图片
               :param file: 附件路径，如果不在当前目录下，要写绝对路径，默认没有附件
               :param ssl: 是否安全链接，默认为安全链接
               :param email_host: smtp服务器地址，默认为qq服务器
               :param ssl_port: 安全链接端口，默认为465
        '''
 
        self.username = "123456789@qq.com"  # 用户名（发件人的邮箱账号或用户名）
        self.passwd = "abcdefghijklmn"  # 16位的SMTP授权码（不含空格）
        self.recv = ['123456789@qq.com'] # 收件人，多个要传list ['a@qq.com','b@qq.com]
        self.title = title  # 邮件标题
        self.content = content  # 邮件正文
        self.file = file  # 绝对路径
        self.image = image  # 图片路径（绝对路径）
        self.email_host = 'smtp.qq.com'  # smtp服务器地址,默认为qq邮箱的服务器
        self.ssl = True,  # 是否安全链接
        self.ssl_port = ssl_port  # 安全链接端口
        # 构造一个MIMEMultipart对象代表邮件本身
        self.message = MIMEMultipart()
 
    # 添加文件到附件
        if self.file:
            file_name = os.path.split(self.file)[-1]  # 只取文件名，不取路径
            try:
                f = open(self.file, 'rb').read()
            except Exception as e:
                traceback.print_exc()
            else:
                att = MIMEText(f, "base64", "utf-8")
                att["Content-Type"] = 'application/octet-stream'
                # base64.b64encode(file_name.encode()).decode()
                new_file_name = '=?utf-8?b?' + base64.b64encode(file_name.encode()).decode() + '?='
                # 处理文件名为中文名的文件
                att["Content-Disposition"] = 'attachment; filename="%s"' % (new_file_name)
                self.message.attach(att)
 
    # 添加图片到附件
        if self.image:
            image_name = os.path.split(self.image)[-1]  # 只取文件名和后缀，不取路径
            try:
                with open(self.image, 'rb') as f:
                    # 图片添加到附件
                    mime = MIMEBase('image', 'image', filename=image_name)
                    mime.add_header('Content-Disposition', 'attachment', filename=image_name)
                    mime.set_payload(f.read())
                    encoders.encode_base64(mime)
                    self.message.attach(mime)
                with open(self.image,'rb') as f:  # 以二进制的方式进行图片读写
                    # 将图片显示在邮件正文中
                    msgimage = MIMEImage(f.read())
                    msgimage.add_header('Content-ID', '<image1>')# 指定文件的Content-ID,<img>,在HTML中图片src将用到
                    self.message.attach(msgimage)
            except Exception as e:
                traceback.print_exc()
 
 
    def send_mail(self):
        self.message.attach(MIMEText(self.content, 'html', 'utf-8'))  # 正文内容   plain代表纯文本,html代表支持html文本
        self.message['From'] = self.username  # 发件人邮箱
        self.message['To'] = ','.join(self.recv)  # 收件人邮箱
        self.message['Subject'] = self.title  # 邮件标题
        self.smtp = smtplib.SMTP_SSL(self.email_host, port=self.ssl_port)
        # 发送邮件服务器的对象
        self.smtp.login(self.username, self.passwd)
        try:
            self.smtp.sendmail(self.username, self.recv, self.message.as_string())
            pass
        except Exception as e:
            resultCode = 0
            traceback.print_exc()
        else:
            result = "邮件发送成功！"
            print(result)
            resultCode = 1
        self.smtp.quit()
        return resultCode  # 定义一个邮件发送结果参数，1为发送成功，0为发送失败。
 
 
if __name__ == "__main__":
    content = '''
               <h1>'''+ sys.argv[1] + '''</h1>
               <p>图片展示：</p>
               <p><img src="cid:image1"></p>
             ''' #通过执行时候的第一个参数输入来确定邮件的内容
    file = os.getcwd() + "/" + "loop.sh" # 附加的文本文件
    image = os.getcwd() +  "/" + "sushe2.png"  # 附件中的图片文件
    m = SendMail("UNTITLE",content=content,file = file,image=image)
    m.send_mail()
    # print(file)
```
执行命令
```shell
python3 filename.py "test"
```
在`Python`中可以利用模块`os`实现和Linux下很多相同功能的操作,`os.getcwd()`就可以获取当前文件的路径.

最后的结果如下图所示,可以将图片修改成自己计算的结果,这样就可以直接通过邮件查看结果,不用再登录服务器去查看结果.


![png](/assets/images/shell/mail-pic.png)


